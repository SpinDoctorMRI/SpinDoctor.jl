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
