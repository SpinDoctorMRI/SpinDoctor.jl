"""
    assemble_flux_matrix(facets, nodes[, weights])

Assemble 3D flux matrix using P1 finite elements.

This function is based on the Matlab function `flux_matrixP1_3D.m` from the Matlab
nodal matrix assembly toolbox by Jan Valdman and Talal Rahman.

https://uk.mathworks.com/matlabcentral/fileexchange/27826-fast-fem-assembly-nodal-elements

Talal Rahman and Jan Valdman: Fast MATLAB assembly of FEM matrices in 2D and 3D: nodal
elements, Applied Mathematics and Computation 219, 7151–7158 (2013).
"""
function assemble_flux_matrix(facets, nodes, weights = 1)
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

    # Enforce symmetry
    sym(x) = (x + x') / 2

    sym(sparse(x[:], y[:], z[:], nnode, nnode))
end

"""
    compute_areas(facets, nodes)

Compute facet areas.
"""
function compute_areas(facets, nodes)
    nfacet = size(facets, 1)

    v1 = nodes[facets[:, 2], :] - nodes[facets[:, 1], :]
    v2 = nodes[facets[:, 3], :] - nodes[facets[:, 1], :]

    matrix_3D = zeros(2, 2, size(facets, 1))
    matrix_3D[1, 1, :] = sum(v1 .* v1, dims = 2)
    matrix_3D[1, 2, :] = sum(v1 .* v2, dims = 2)
    matrix_3D[2, 1, :] = sum(v1 .* v2, dims = 2)
    matrix_3D[2, 2, :] = sum(v2 .* v2, dims = 2)

    [√abs(det(matrix_3D[:, :, i])) / 2 for i = 1:nfacet]
end
