function compute_int_J(params, α)
    (; r, D, W, dim) = params
    m = length(D)
    R = r[m]

    λ = α^2

    if α > 1e-6
        V, _ = compute_v(params, α, 0, R)
        J = W[m] * V / λ / R
    else
        J = 1 / dim
    end

    J
end
