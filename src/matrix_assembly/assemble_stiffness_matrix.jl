"""
    assemble_stiffness_matrix(elements, nodes[, D])

Assemble 3D stiffness matrix using P1 finite elements.
A diffusion tensor `D` may be provided.

This function is based on the Matlab function `stiffness_matrixP1_3D.m` from the Matlab
nodal matrix assembly toolbox by Jan Valdman and Talal Rahman.

https://uk.mathworks.com/matlabcentral/fileexchange/27826-fast-fem-assembly-nodal-elements

Talal Rahman and Jan Valdman: Fast MATLAB assembly of FEM matrices in 2D and 3D: nodal
elements, Applied Mathematics and Computation 219, 7151–7158 (2013).
"""
function assemble_stiffness_matrix(elements, nodes, D = I(3))
    nelement = size(elements, 1)
    nodes = permutedims(nodes[elements, :], [3, 2, 1])

    integration_point = [1 // 4; 1 // 4; 1 // 4]
    ∇φ, detjac, _ = compute_∇φ(nodes, integration_point)

    volumes = abs.(detjac[:]) / 6
    ∇φ = ∇φ[:, :, 1, :]

    z = zeros(4, 4, nelement)
    for ielement = 1:nelement
        z[:, :, ielement] = volumes[ielement] * ∇φ[:, :, ielement]' * D * ∇φ[:, :, ielement]
    end
    y = reshape(repeat(elements, 1, 4)', 4, 4, nelement)
    x = permutedims(y, [2, 1, 3])

    # Enforce symmetry
    sym(x) = (x + x') / 2

    sym(sparse(x[:], y[:], z[:]))
end

"""
    compute_∇φ(nodes, points)

Compute the gradients of the finite element basis functions φ.
"""
function compute_∇φ(nodes, points)
    nelement = size(nodes, 3)
    npoint = size(points, 2)
    dshape = shapeder(points)
    ∇φ = zeros(3, 4, npoint, nelement)
    detjac = zeros(npoint, nelement)
    jac = zeros(3, 3, npoint, nelement)
    for ipoint = 1:npoint
        for ielement = 1:nelement
            j = dshape[:, :, ipoint] * nodes[:, :, ielement]'
            jinv = inv(j)
            detjac[ipoint, ielement] = abs(det(j))
            jac[:, :, ipoint, ielement] = j
            ∇φ[:, :, ipoint, ielement] = jinv * dshape[:, :, ipoint]
        end
    end
    ∇φ, detjac, jac
end

"""
    shapeder(points)

Derivative of shape functions with respect to reference coordinates (`ξ₁`, `ξ₂`, `ξ₃`).
"""
function shapeder(points)
    npoint = size(points, 2)
    dshape = [
        1.0 0.0 0.0 -1.0
        0.0 1.0 0.0 -1.0
        0.0 0.0 1.0 -1.0
    ]
    dshape = repeat(dshape, 1, 1, npoint)
end
