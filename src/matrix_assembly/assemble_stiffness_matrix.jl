""" Assemble 3D stiffness matrix using P1 finite elements. """
function assemble_stiffness_matrix(elements, nodes, weights=1)
    nelement = size(elements, 1)
    nnode = size(nodes, 1)
    @assert length(size(weights)) ≤ 1
    nweight = size(weights, 1)
    nodes = permutedims(nodes[elements, :], [3, 2, 1])

    integration_point = [1//4; 1//4; 1//4]
    ∇φ, detjac, _ = compute_∇φ(nodes, integration_point)

    volumes = abs.(detjac[:]) / 6
    ∇φ = ∇φ[:, :, 1, :]

    z = zeros(4, 4, nelement)
    for ielement = 1:nelement
        z[:, :, ielement] = volumes[ielement] * weights * ∇φ[:, :, ielement]' * ∇φ[:, :, ielement]
    end
    y = reshape(repeat(elements, 1, 4)', 4, 4, nelement)
    x = permutedims(y, [2, 1, 3])

    # Assure symmetry
    sym(x) = (x + x') / 2

    sym(sparse(x[:], y[:], z[:]))
end

""" Gradients of the finite element basis functions φ. """
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

""" Derivative of shape functions with respect to reference coordinates (ξ1, ξ2, ξ3). """
function shapeder(points)
    npoint = size(points, 2)
    dshape = [
        1. 0. 0. -1.
        0. 1. 0. -1.
        0. 0. 1. -1.
    ]
    dshape = repeat(dshape, 1, 1, npoint)
end
