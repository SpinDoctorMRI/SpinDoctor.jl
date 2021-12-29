function compute_int_K(params, α₁, n₁, α₂, n₂)
    (; r, D, dim) = params

    m = length(D)
    R = r[m]

    dr = r[m] * 1e-12
    r₁ = [0; r[1:m-1]] .+ dr
    r₂ = r[1:m] .- dr

    K = 0
    if abs(n₁ - n₂) == 1
        λ₁ = α₁^2
        λ₂ = α₂^2
        ν₁ = compute_ν(n₁, dim)
        ν₂ = compute_ν(n₂, dim)

        # Computation at point r₁
        V₁, dV₁ = compute_v(params, α₁, n₁, r₁)
        V₂, dV₂ = compute_v(params, α₂, n₂, r₁)
        I₁a =
            (
                (λ₁ + λ₂) * D .* r₁ .^ (dim - 1) -
                (ν₁ + ν₂ - (dim - 1)) * D .^ 2 .* r₁ .^ (dim - 3)
            ) .* V₁ .* V₂
        I₁b =
            (
                (λ₂ - λ₁) * D .* r₁ .^ dim -
                (ν₁ - ν₂ - (dim - 1)) * D .^ 2 .* r₁ .^ (dim - 2)
            ) .* dV₁ .* V₂
        I₁c =
            (
                (λ₁ - λ₂) * D .* r₁ .^ dim -
                (ν₂ - ν₁ - (dim - 1)) * D .^ 2 .* r₁ .^ (dim - 2)
            ) .* V₁ .* dV₂
        I₁d = 2 * D .^ 2 .* r₁ .^ (dim - 1) .* dV₁ .* dV₂
        I₁ = I₁a + I₁b + I₁c + I₁d

        # Computation at point r₂
        V₁, dV₁ = compute_v(params, α₁, n₁, r₂)
        V₂, dV₂ = compute_v(params, α₂, n₂, r₂)
        I₂a =
            (
                (λ₁ + λ₂) * D .* r₂ .^ (dim - 1) -
                (ν₁ + ν₂ - (dim - 1)) * D .^ 2 .* r₂ .^ (dim - 3)
            ) .* V₁ .* V₂
        I₂b =
            (
                (λ₂ - λ₁) * D .* r₂ .^ dim -
                (ν₁ - ν₂ - (dim - 1)) * D .^ 2 .* r₂ .^ (dim - 2)
            ) .* dV₁ .* V₂
        I₂c =
            (
                (λ₁ - λ₂) * D .* r₂ .^ dim -
                (ν₂ - ν₁ - (dim - 1)) * D .^ 2 .* r₂ .^ (dim - 2)
            ) .* V₁ .* dV₂
        I₂d = 2 * D .^ 2 .* r₂ .^ (dim - 1) .* dV₁ .* dV₂
        I₂ = I₂a + I₂b + I₂c + I₂d
    end

    K = sum(I₂ - I₁) / (λ₁ - λ₂)^2 / R^(dim + 1)
end
