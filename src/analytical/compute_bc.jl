function compute_bc(params, α, n)
    (; r, D, W, dim) = params
    m = length(D)

    bc = zeros(m, 2)
    bc[:, 1] .= 1

    if α == 0
        return bc
    end

    ρ = r[1:m-1] ./ .√D[1:m-1]
    ρt = r[1:m-1] ./ .√D[2:m]

    for i = 1:m-1
        J, Y, dJ, dY = compute_JY(α * ρ[i], n, dim)
        J1, Y1, dJ1, dY1 = compute_JY(α * ρt[i], n, dim)

        A11 = -√(D[i] / D[i+1]) * dJ * Y1 + J * dY1 + √D[i] / W[i] * α * dJ * dY1
        A12 = -√(D[i] / D[i+1]) * dY * Y1 + Y * dY1 + √D[i] / W[i] * α * dY * dY1
        A21 = √(D[i] / D[i+1]) * dJ * J1 - J * dJ1 - √D[i] / W[i] * α * dJ * dJ1
        A22 = √(D[i] / D[i+1]) * dY * J1 - Y * dJ1 - √D[i] / W[i] * α * dY * dJ1

        z = α * ρt[i]
        if dim == 2
            Q = z * π / 2
        else
            Q = z^2
        end
        A = Q * [A11 A12; A21 A22]
        bc[i+1, :] = A * bc[i, :]
    end

    bc
end
