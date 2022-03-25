"""
    assemble_mass_matrix(elements, volumes[, weights])

Assemble 3D mass matrix using P1 finite elements.

This function is based on the Matlab function `mass_matrixP1_3D.m` from the Matlab
nodal matrix assembly toolbox by Jan Valdman and Talal Rahmnan:
    
https://uk.mathworks.com/matlabcentral/fileexchange/27826-fast-fem-assembly-nodal-elements
    
Talal Rahman and Jan Valdman: Fast MATLAB assembly of FEM matrices in 2D and 3D: nodal
elements, Applied Mathematics and Computation 219, 7151–7158 (2013).
"""
function assemble_mass_matrix(elements, volumes, weights = 1)
    x = kron(ones(1, 4), elements)
    y = kron(elements, ones(1, 4))
    if length(weights) == 1 || length(weights) == size(elements, 1)
        # P0 weights
        z = kron(weights .* volumes, reshape((ones(4, 4) + I(4)) / 20, 1, 16))
    else
        # P1 weights
        M1 = [6 2 2 2; 2 2 1 1; 2 1 2 1; 2 1 1 2] / 120
        M2 = M1[[4, 1, 2, 3], [4, 1, 2, 3]]
        M3 = M2[[4, 1, 2, 3], [4, 1, 2, 3]]
        M4 = M3[[4, 1, 2, 3], [4, 1, 2, 3]]
        z = (
            kron(volumes .* weights[elements[:, 1]], reshape(M1, 1, 16)) +
            kron(volumes .* weights[elements[:, 2]], reshape(M2, 1, 16)) +
            kron(volumes .* weights[elements[:, 3]], reshape(M3, 1, 16)) +
            kron(volumes .* weights[elements[:, 4]], reshape(M4, 1, 16))
        )
    end

    # Enforce symmetry
    sym(x) = (x + x') / 2

    sym(sparse(x[:], y[:], z[:]))
end
