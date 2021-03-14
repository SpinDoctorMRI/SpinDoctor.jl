""" Assemble 3D flux matrix using P1 finite elements. """
function assemble_flux_matrix(facets, nodes, weights = 1)
    nfacet = size(facets, 1)
    nnode = size(nodes, 1)
    areas = compute_areas(facets, nodes)

    x = kron(ones(1, 3), facets)
    y = kron(facets, ones(1, 3))
    if length(weights) == 1 || length(weights) == size(facets, 1)
        # P0 weights
        z = kron(weights .* areas, reshape((ones(3, 3) + I(3)) / 12, 1, 9))
    else
        # P1 weights
        M1 = [6 2 2; 2 2 1; 2 1 2] / 60
        M2 = M1[[3, 1, 2], [3, 1, 2]]
        M3 = M2[[3, 1, 2], [3, 1, 2]]

        z = (
            kron(areas .* weights[facets[:, 1]], reshape(M1, 1, 9)) +
            kron(areas .* weights[facets[:, 2]], reshape(M2, 1, 9)) +
            kron(areas .* weights[facets[:, 3]], reshape(M3, 1, 9))
        )
    end

    # Assure symmetry
    sym(x) = (x + x') / 2

    sym(sparse(x[:], y[:], z[:], nnode, nnode))
end

""" Compute facet areas. """
function compute_areas(facets, nodes)
    nfacet = size(facets, 1)

    v1 = nodes[facets[:, 2], :] - nodes[facets[:, 1], :]
    v2 = nodes[facets[:, 3], :] - nodes[facets[:, 1], :]

    matrix_3D = zeros(2, 2, size(facets, 1))
    matrix_3D[1, 1, :] = sum(v1 .* v1, dims = 2)
    matrix_3D[1, 2, :] = sum(v1 .* v2, dims = 2)
    matrix_3D[2, 1, :] = sum(v1 .* v2, dims = 2)
    matrix_3D[2, 2, :] = sum(v2 .* v2, dims = 2)

    areas = ([√abs(det(matrix_3D[:, :, i])) / 2 for i = 1:nfacet])
end
