"""
    get_mesh_surface(points, facets)

Compute surface areas, centers and normals for each facet.
"""
function get_mesh_surface(points, facets)
    dim = size(points, 1)
    nfacet = size(facets, 2)

    # Facets
    x = points[:, facets]

    # Facet centers
    centers = reshape(mean(x; dims = 2), dim, nfacet)

    # Facet normals
    if dim == 2
        normals = reshape(
            mapslices(x; dims = [1, 2]) do x
                v = x[:, 1] - x[:, 2]
                [v[2], -v[1]]
            end,
            2,
            nfacet,
        )
    elseif dim == 3
        normals = reshape(
            mapslices(x -> 1 / 2 * (x[:, 1] - x[:, 2]) × (x[:, 3] - x[:, 2]), x; dims = [1, 2]),
            3,
            nfacet,
        )
    end

    # Facet areas
    areas = reshape(mapslices(norm, normals; dims = 1), :)

    # Normalize normals
    normals = normals ./ sqrt.(sum(abs2, normals; dims = 1))

    areas, centers, normals
end
