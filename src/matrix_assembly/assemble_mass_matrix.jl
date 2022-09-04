"""
    assemble_mass_matrix(elements, nodes[, weights])

Assemble 3D mass matrix using P1 finite elements.

This function is based on the Matlab function `mass_matrixP1_3D.m` from the Matlab
nodal matrix assembly toolbox by Jan Valdman and Talal Rahmnan:

<https://uk.mathworks.com/matlabcentral/fileexchange/27826-fast-fem-assembly-nodal-elements>

[Rahman2013](@cite)
"""
function assemble_mass_matrix(elements, nodes, weights = 1)
    nnode = size(nodes, 2)
    nweight = length(weights)
    nbasis, nelement = size(elements)
    if nbasis == 2
        volumes = compute_lengths(elements, nodes)
    elseif nbasis == 3
        volumes = compute_areas(elements, nodes)
    elseif nbasis == 4
        volumes, = get_mesh_volumes(nodes, elements)
    end
    x = kron(ones(nbasis), elements)'
    y = kron(elements, ones(nbasis))'
    if nweight == 1 || nweight == nelement
        # P0 weights
        if nbasis == 2
            b = (ones(2, 2) + I(2)) / 3
        elseif nbasis == 3
            b = (ones(3, 3) + I(3)) / 12
        elseif nbasis == 4
            b = (ones(4, 4) + I(4)) / 20
        end
        z = kron(weights .* volumes, reshape(b, 1, :))
    else
        # P1 weights
        if nbasis == 2
        elseif nbasis == 3
            M1 = [6 2 2; 2 2 1; 2 1 2] / 60
        elseif nbasis == 4
            M1 = [6 2 2 2; 2 2 1 1; 2 1 2 1; 2 1 1 2] / 120
        end
        z = sum(1:nbasis) do i
            kron(
                volumes .* weights[elements[i, :]],
                reshape(circshift(M1, (i - 1, i - 1)), 1, :),
            )
        end
    end

    # Enforce symmetry
    sym(x) = (x + x') / 2

    sym(sparse(x[:], y[:], z[:], nnode, nnode))
end
