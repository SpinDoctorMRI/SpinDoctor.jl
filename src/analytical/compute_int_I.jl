function compute_int_I(params, α₁, α₂, n)
    @unpack r, D, dim = params
    m = length(D)
    R = r[m]

    dr = r[m] * 1e-12
    r₂ = r[1:m] .- dr
    r₁ = [0; r[1:m-1]] .+ dr

    # Check if α₁ is not equal to α₂
    if abs(α₁ - α₂) > 1e-8
        λ₁ = α₁^2
        λ₂ = α₂^2

        # Computation at point r₁
        V₁, dV₁ = compute_v(params, α₁, n, r₁)
        V₂, dV₂ = compute_v(params, α₂, n, r₁)
        I₁ = D .* r₁ .^ (dim - 1) .* (dV₁ .* V₂ - dV₂ .* V₁)

        # Computation at point r₂
        V₁, dV₁ = compute_v(params, α₁, n, r₂)
        V₂, dV₂ = compute_v(params, α₂, n, r₂)
        I₂ = D .* r₂ .^ (dim - 1) .* (dV₁ .* V₂ - dV₂ .* V₁)

        I = (I₂ - I₁) / (λ₁ - λ₂) / R^dim
    else
        α = α₁
        λ = α^2
        nu = compute_ν(n, dim)

        if α > 1e-6
            V, dV = compute_v(params, α, n, r₁)
            I₁ = (
                D .* r₁ .^ dim .* dV .^ 2 +
                (λ * r₁ .^ dim - D .* r₁ .^ (dim - 2) * nu) .* V .^ 2 +
                (dim - 2) * D .* r₁ .^ (dim - 1) .* dV .* V
            )

            V, dV = compute_v(params, α, n, r₂)
            I₂ = (
                D .* r₂ .^ dim .* dV .^ 2 +
                (λ * r₂ .^ dim - D .* r₂ .^ (dim - 2) * nu) .* V .^ 2 +
                (dim - 2) * D .* r₂ .^ (dim - 1) .* dV .* V
            )

            I = (I₂ - I₁) / (2 * λ * R^dim)
        else
            I = (r₂ .^ dim - r₁ .^ dim) / R^dim / dim
        end
    end

    I
end