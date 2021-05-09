function compute_v(params, α, n, x, bc = compute_bc(params, α, n))
    @unpack r, D, dim = params

    x = [x;]

    m = length(D)
    v = ones(size(x))
    dv = zeros(size(x))
    if α > 0
        r₁ = [0; r[1:m]]
        r₁[1] -= 1e-8
        r₁[m+1] += 1e-8
        for i = 1:m
            inds = r₁[i] .< x .≤ r₁[i+1]
            J, Y, dJ, dY = compute_JY(α * x[inds] / √D[i], n, dim)
            if abs(bc[i, 2]) > 1e-12
                v[inds] = bc[i, 1] * J + bc[i, 2] * Y
                dv[inds] = α / √D[i] * (bc[i, 1] * dJ + bc[i, 2] * dY)
            else
                v[inds] = bc[i, 1] * J
                dv[inds] = α / √D[i] * bc[i, 1] * dJ
            end
        end
    end

    v, dv
end