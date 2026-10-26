struct FEMesh{T,dim}
    point_map::Vector{Vector{Int}}
    points::Vector{Matrix{T}}
    facets::Matrix{Matrix{Int}}
    elements::Vector{Matrix{Int}}
end

function FEMesh(; point_map, points, facets, elements)
    T = eltype(first(points))
    dim = size(first(points), 1)
    FEMesh{T,dim}(point_map, points, facets, elements)
end

"""
    split_mesh(mesh_all)

Split global mesh into compartments.
"""
function split_mesh(mesh_all)
    # Extract mesh_all in global numbering system (tag: "all")
    points_all = mesh_all.points
    facets_all = mesh_all.facets
    facetmarkers_all = mesh_all.facetmarkers
    elements_all = mesh_all.elements
    elementmarkers_all = mesh_all.elementmarkers

    dim, npoint = size(mesh_all.points)
    ncompartment = maximum(elementmarkers_all)
    nboundary = maximum(facetmarkers_all)
    elements = [elements_all[:, elementmarkers_all .== icmpt] for icmpt = 1:ncompartment]

    # Map global to local points
    point_map = sort.(unique.(elements))
    points = [points_all[:, point_map[icmpt]] for icmpt = 1:ncompartment]

    # Assign each boundary facet to the compartments that have it as an element face
    # (an interface belongs to both of its compartments), and renumber locally
    facets = [zeros(Int, dim, 0) for _ = 1:ncompartment, _ = 1:nboundary]
    point_map_inv = zeros(Int, npoint)
    for icmpt = 1:ncompartment
        point_map_inv[point_map[icmpt]] = 1:length(point_map[icmpt])
        inside = facet_elements(elements[icmpt], facets_all) .> 0
        for iboundary = 1:nboundary
            sel = inside .& (facetmarkers_all .== iboundary)
            if any(sel)
                facets[icmpt, iboundary] = point_map_inv[facets_all[:, sel]]
            end
        end
        elements[icmpt] = point_map_inv[elements[icmpt]]
        point_map_inv[point_map[icmpt]] .= 0
    end

    FEMesh(; point_map, points, facets, elements)
end

"""
    split_field(mesh, ξ)

Split the field `ξ` into a vector containing a `view` of each compartment. 
"""
function split_field(mesh, ξ)
    ncompartment = length(mesh.points)
    npoint_cmpts = size.(mesh.points, 2)
    inds = [0; cumsum(npoint_cmpts[1:end])]
    [@view(ξ[(1+inds[i]):inds[i+1]]) for i = 1:ncompartment]
end

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

    onplane(f, d, value) =
        all(a -> abs(points[d, facets[a, f]] - value) ≤ atol, axes(facets, 1))
    face(d, value) = facets[:, [f for f ∈ axes(facets, 2) if onplane(f, d, value)]]

    names = (:xmin, :xmax, :ymin, :ymax, :zmin, :zmax)[1:2dim]
    faces = [face(d, side == 1 ? pmin[d] : pmax[d]) for d = 1:dim for side = 1:2]
    NamedTuple{names}(faces)
end

"""
    facet_elements(elements, facets)

For each facet (column of `facets`), find an element (column of `elements`) that has the
facet as one of its faces. The result is zero for facets that are not a face of any element.
"""
facet_elements(elements, facets) = facet_elements(elements, facets, Val(size(facets, 1)))

function facet_elements(elements, facets, ::Val{d}) where {d}
    nbasis = size(elements, 1)
    face(f, j) = sort(SVector{d}(ntuple(k -> f[k, j], Val(d))))

    # Look up the faces of all elements (drop one vertex at a time) among the facets
    lookup = Dict(face(facets, j) => j for j ∈ axes(facets, 2))
    ielements = zeros(Int, size(facets, 2))
    for e ∈ axes(elements, 2), skip = 1:nbasis
        f = sort(SVector{d}(ntuple(k -> elements[k < skip ? k : k + 1, e], Val(d))))
        j = get(lookup, f, 0)
        j > 0 && (ielements[j] = e)
    end

    # Repeated facets share the element of the one in the lookup table
    for j ∈ axes(facets, 2)
        ielements[j] = ielements[lookup[face(facets, j)]]
    end

    ielements
end

"""
    get_cmpt_volumes(mesh)

Get compartment volumes.
"""
function get_cmpt_volumes(mesh)
    map((p, e) -> sum(get_mesh_volumes(p, e)[1]), mesh.points, mesh.elements)
end

"""
    get_mesh_volumes(points, elements)

Get the volumes and centers of the elements of a mesh.
"""
function get_mesh_volumes(points, elements)
    dim = size(points, 1)
    n = size(elements, 1)
    volumes = _element_measures(points, elements, Val(dim), Val(n))
    centers = reshape(mean(reshape(points[:, elements], dim, n, :); dims = 2), dim, :)
    volumes, centers
end

_element_measures(points, elements, dim, n) =
    [simplex_measure(simplex(points, elements, e, dim, n)) for e ∈ axes(elements, 2)]

"""
    get_mesh_surface(points, facets)

Compute surface areas, centers and unit normals for each facet. The orientation of the
normals is given by the order of the facet points (see
[`get_mesh_surfacenormals`](@ref) for outward normals).
"""
get_mesh_surface(points, facets) = _get_mesh_surface(points, facets, Val(size(points, 1)))

function _get_mesh_surface(points::AbstractMatrix{T}, facets, ::Val{dim}) where {T,dim}
    nfacet = size(facets, 2)
    areas = zeros(T, nfacet)
    centers = zeros(T, dim, nfacet)
    normals = zeros(T, dim, nfacet)
    for f = 1:nfacet
        x = simplex(points, facets, f, Val(dim), Val(dim))
        if dim == 2
            v = x[:, 1] - x[:, 2]
            n = SVector(v[2], -v[1])
        else
            n = (x[:, 1] - x[:, 2]) × (x[:, 3] - x[:, 2]) / 2
        end
        areas[f] = norm(n)
        normals[:, f] = n / areas[f]
        centers[:, f] = sum(x; dims = 2) / dim
    end
    areas, centers, normals
end

"""
    get_mesh_surfacenormals(points, elements, facets)

Compute areas, centers and outward unit normals of the boundary `facets` of a mesh
compartment with `points` and `elements`.
"""
function get_mesh_surfacenormals(points, elements, facets)
    dim = size(points, 1)
    areas, facet_centers, normals = get_mesh_surface(points, facets)

    # The element of each facet lies on the inner side of the facet
    ielements = facet_elements(elements, facets)
    all(>(0), ielements) || error("Some facets are not faces of the elements")
    element_centers = reshape(mean(points[:, elements[:, ielements]]; dims = 2), dim, :)

    # Orient the normals from the element center towards the facet center
    orientations = sign.(sum((facet_centers - element_centers) .* normals; dims = 1))
    normals = orientations .* normals

    areas, facet_centers, normals
end
