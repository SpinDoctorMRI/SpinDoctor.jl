"""
    assemble_stiffness_matrix(elements, nodes[, D])

Assemble 3D stiffness matrix using P1 finite elements.
A diffusion tensor `D` may be provided.

This function is based on the Matlab function `stiffness_matrixP1_3D.m` from the Matlab
nodal matrix assembly toolbox by Jan Valdman and Talal Rahman.

<https://uk.mathworks.com/matlabcentral/fileexchange/27826-fast-fem-assembly-nodal-elements>

[Rahman2013](@cite)
"""
function assemble_stiffness_matrix(elements, nodes, D = I(3))
    dim, nnode = size(nodes)
    nbasis, nelement = size(elements)

    if dim == 2
        integration_point = [1 // 3; 1 // 3]
    elseif dim == 3
        integration_point = [1 // 4; 1 // 4; 1 // 4]
    end

    ∇φ, detjac, _ = compute_∇φ(nodes[:, elements], integration_point)
    volumes = abs.(detjac[:]) / factorial(dim)
    ∇φ = ∇φ[:, :, 1, :]

    z = zeros(nbasis, nbasis, nelement)
    for ielement = 1:nelement
        z[:, :, ielement] = volumes[ielement] * ∇φ[:, :, ielement]' * D * ∇φ[:, :, ielement]
    end
    y = reshape(repeat(elements, nbasis), nbasis, nbasis, nelement)
    x = permutedims(y, [2, 1, 3])

    # Enforce symmetry
    sym(x) = (x + x') / 2

    sym(sparse(x[:], y[:], z[:], nnode, nnode))
end

"""
    compute_∇φ(nodes, points)

Compute the gradients of the finite element basis functions φ.
"""
function compute_∇φ(nodes, points)
    dim, nbasis, nelement = size(nodes)
    npoint = size(points, 2)
    dshape = shapeder(points)
    ∇φ = zeros(dim, nbasis, npoint, nelement)
    detjac = zeros(npoint, nelement)
    jac = zeros(dim, dim, npoint, nelement)
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
    dim = size(points, 1)
    npoint = size(points, 2)
    dshape = [I(dim) -ones(dim)]
    repeat(dshape, 1, 1, npoint)
end
