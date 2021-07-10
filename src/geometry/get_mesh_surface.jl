"""
    get_mesh_surface(points, facets)
    
Compute surface areas, centers and normals for each facet.
"""
function get_mesh_surface(points, facets)

    # Number of facets
    nfacet = size(facets, 2)

    # Facets
    x = points[:, facets]

    # Facet centers
    centers = reshape(mean(x, dims = 2), 3, nfacet)

    # Facet normals
    normals = reshape(
        mapslices(x -> (x[:, 1] - x[:, 2]) × (x[:, 3] - x[:, 2]), x, dims = [1, 2]),
        3,
        nfacet,
    )

    # Facet areas
    areas = 1 / 2 * reshape(mapslices(norm, normals, dims = 1), nfacet)
    total_area = sum(areas)

    # Normalize normals
    normals = normals ./ 2areas'

    total_area, areas, centers, normals
end
