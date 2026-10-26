"""
    compartment_mesh(mesh, icmpt)

Extract compartment `icmpt` of `mesh` as a mesh of its own, with one compartment. The
boundary numbering is kept (`facets` has one row and the same columns as `mesh.facets`), and
`point_map` still maps to the points of the global mesh.

A model of the compartment alone can then be built with coefficients for one compartment,
for example to solve a problem in the extracellular space only.
"""
function compartment_mesh(mesh::FEMesh, icmpt)
    FEMesh(;
        point_map = [mesh.point_map[icmpt]],
        points = [mesh.points[icmpt]],
        facets = mesh.facets[icmpt:icmpt, :],
        elements = [mesh.elements[icmpt]],
    )
end

"""
    box_faces(mesh, icmpt = length(mesh.points); atol = 1e-8 * diameter)

Find the boundary facets of compartment `icmpt` (by default the last one, which is the ECS
in the setups with an ECS) lying on the faces of its bounding box. This applies to
axis-aligned box-shaped compartments, such as a [`BoxECS`](@ref) that has not been bent or
twisted.

Returns a named tuple with the facets (in the local numbering of the compartment) on each
face: `(; xmin, xmax, ymin, ymax)` in 2D and `(; xmin, xmax, ymin, ymax, zmin, zmax)` in 3D.
The points on a face are given by `unique(facets)`.
"""
function box_faces(
    mesh::FEMesh{T,dim},
    icmpt = length(mesh.points);
    atol = nothing,
) where {T,dim}
    points = mesh.points[icmpt]
    facets = reduce(hcat, mesh.facets[icmpt, :])
    pmin = minimum(points; dims = 2)
    pmax = maximum(points; dims = 2)
    isnothing(atol) && (atol = 1e-8 * norm(pmax - pmin))

    onplane(f, d, value) = all(a -> abs(points[d, facets[a, f]] - value) ≤ atol, axes(facets, 1))
    face(d, value) = facets[:, [f for f ∈ axes(facets, 2) if onplane(f, d, value)]]

    names = (:xmin, :xmax, :ymin, :ymax, :zmin, :zmax)[1:2dim]
    faces = [face(d, side == 1 ? pmin[d] : pmax[d]) for d = 1:dim for side = 1:2]
    NamedTuple{names}(faces)
end
