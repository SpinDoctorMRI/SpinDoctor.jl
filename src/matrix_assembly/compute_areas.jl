"""
    compute_areas(facets, nodes)

Compute facet areas.
"""
function compute_areas(facets, nodes)
    nfacet = size(facets, 2)

    v1 = nodes[:, facets[2, :]] - nodes[:, facets[1, :]]
    v2 = nodes[:, facets[3, :]] - nodes[:, facets[1, :]]

    matrix_3D = zeros(2, 2, nfacet)
    matrix_3D[1, 1, :] = sum(v1 .* v1; dims = 1)
    matrix_3D[1, 2, :] = sum(v1 .* v2; dims = 1)
    matrix_3D[2, 1, :] = sum(v1 .* v2; dims = 1)
    matrix_3D[2, 2, :] = sum(v2 .* v2; dims = 1)

    [sqrt(abs(det(matrix_3D[:, :, i]))) / 2 for i = 1:nfacet]
end
