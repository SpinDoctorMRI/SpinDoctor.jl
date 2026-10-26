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
