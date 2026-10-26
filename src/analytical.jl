# This module ("Analytical") is based on the following articles and corresponding MATLAB code:
#   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
#       Rev. Mod.Phys. 79, 1077 (2007)
#   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted
#       diffusion inmultilayered structures,
#       J. Magn. Reson. 205, 181-195 (2010).

"""
    solve_analytical_laplace(; ρ, r, D, W, T₂, γ, dim, eiglim, eigstep)

Solve the analytical radial Laplace eigenvalue problem of multilayered plates (`dim = 1`),
cylinders (`dim = 2`) or spheres (`dim = 3`), with the coefficients of
[`analytical_coefficients`](@ref), for eigenvalues up to `eiglim`, searched with steps
`eigstep`. Return the coefficients together with the eigendecomposition.
"""
function solve_analytical_laplace(; ρ, r, D, W, T₂, γ, dim, eiglim, eigstep)
    laplace = (; ρ, r, D, W, T₂, γ, dim, eiglim, eigstep)

    # TODO: Add support for dim = 1 (PlateSetup)
    L = length(D)

    # Compute radial and angular eigenvalues (λ and ν), with α^2 = λ and n^2 = ν for
    # cylinders and n(n+1) = ν for spheres
    dα = √eigstep * 1e3
    α_min = dα / 100
    α_max = √eiglim * 1e3
    α, n = find_α(laplace, α_min, α_max, dα)

    N = length(α)
    λ = α .^ 2
    Λ = diagm(λ)

    β = zeros(N)
    J = zeros(N)
    for m = 1:N
        β[m] = compute_β(laplace, α[m], n[m])
        if n[m] == 0
            J[m] = compute_int_J(laplace, α[m])[1]
        end
    end

    if dim == 1
        coeff = √2
    elseif dim == 2
        coeff = 2
    elseif dim == 3
        coeff = √6
    end
    U = coeff * J .* β

    B = zeros(N, N)
    for m1 = 1:N
        for m2 = m1:N
            if dim == 1
                K = compute_int_K(laplace, α[m1], n[m1], α[m2], n[m2])
                B[m1, m2] = 2 * β[m1] * β[m2] * K
                B[m2, m1] = B[m1, m2]
            elseif abs(n[m1] - n[m2]) == 1
                K = compute_int_K(laplace, α[m1], n[m1], α[m2], n[m2])
                if dim == 2
                    ε = √(1 + (n[m1] == 0) + (n[m2] == 0))
                elseif dim == 3
                    ε = (n[m1] + n[m2] + 1) / √((2 * n[m1] + 1) * (2 * n[m2] + 1))
                end
                B[m1, m2] = ε * β[m1] * β[m2] * K
                B[m2, m1] = B[m1, m2]
            end
        end
    end

    Bri = zeros(N, N, L)
    for m1 = 1:N
        for m2 = m1:N
            if dim == 1
                int_I = compute_int_I(laplace, α[m1], α[m2], 0)
                Bri[m1, m2, :] = 2 * β[m1] * β[m2] * int_I
                Bri[m2, m1, :] = Bri[m1, m2, :]
            elseif n[m1] == n[m2]
                int_I = compute_int_I(laplace, α[m1], α[m2], n[m1])
                Bri[m1, m2, :] = 2 * β[m1] * β[m2] * int_I
                Bri[m2, m1, :] = Bri[m1, m2, :]
            end
        end
    end

    Br = zeros(N, N)
    for i = 1:L
        Br .+= Bri[:, :, i] / T₂[i]
    end

    merge(laplace, (; Λ, β, U, B, Bri, Br, α, n))
end

"""
    solve_analytical_mf(laplace, volumes, gradient)

Compute the signal in multilayered plates, cylinders or spheres using an analytical matrix
formalism solution, from the eigendecomposition `laplace` of
[`solve_analytical_laplace`](@ref) and the compartment `volumes`. Only implemented for
[`PGSE`](@ref) sequences.

This function is based on the following articles and corresponding code:
    [1] [Grebenkov2007](@cite)
    [2] [Grebenkov2010a](@cite)
"""
function solve_analytical_mf(laplace, volumes, gradient::ScalarGradient)
    (; r, ρ, γ, Λ, Br, B, U) = laplace

    gradient.profile isa PGSE || error("Only implemeted for PGSE")

    w = γ * gradient.amplitude * 1e12 * r[end]
    δ = gradient.profile.δ * 1e-6
    Δ = gradient.profile.Δ * 1e-6

    U₀ = complex(U)
    U = complex(U)

    expmv!(-δ, Λ + Br + im * w * B, U)
    expmv!(-(Δ - δ), Λ + Br, U)
    expmv!(-δ, Λ + Br - im * w * B, U)

    U₀'U * volumes'ρ
end

function find_α(params, α_min, α_max, dα)
    (; dim) = params

    if dim == 1
        α = find_zeros(α -> α_func(params, α, 0), α_min, α_max)
        n = zeros(Int, size(α))
    else
        α = Float64[]
        n = Int[]
        α₀ = [0.0]
        j = 0
        # pl = nothing
        # αα = LinRange(α_min, α_max, 1000)
        while !isempty(α₀) && α₀[1] < α_max
            @info "Computing zeros of Fₙ(α) for n = $j"
            α₀ = find_zeros(α -> α_func(params, α, j), α_min, α_max)
            append!(α, α₀)
            append!(n, fill(j, size(α₀)))
            # if j == 0
            #     pl = lines(αα,  α -> α_func(params, α, j))
            #     display(pl)
            # else
            #     lines!(αα, α -> α_func(params, α, j))
            # end
            j = j + 1
        end

        inds = sortperm(α)
        α = α[inds]
        n = n[inds]
    end

    # Number of layers
    L = length(params.D)

    inds_keep = dα / 100 .< α .< α_max
    α = α[inds_keep]
    n = n[inds_keep]

    if params.W[L] < 1e-12 && (length(params.W) == L || params.W[end] < 1e-12)
        # Add ground eigenmode explicitly
        pushfirst!(α, 0)
        pushfirst!(n, 0)
    end

    α, n
end

function α_func(params, α, n)
    (; r, D, W, dim) = params

    m = length(D)
    bc1 = zeros(size(α))
    bc2 = zeros(size(α))

    for i = 1:length(α)
        bc = compute_bc(params, α[i], n)
        bc1[i] = bc[m, 1]
        bc2[i] = bc[m, 2]
    end

    rD = r[m] / √D[m]

    if dim == 1
        f1 = -sqrt(D[m]) * α .* sin(α * rD) + W[m] * cos(α * rD)
        f2 = sqrt(D[m]) * α .* cos(α * rD) + W[m] * sin(α * rD)
    elseif dim == 2
        f1 =
            √D[m] * α .* (besselj(n - 1, α * rD) - besselj(n + 1, α * rD)) / 2 +
            W[m] * besselj(n, α * rD)
        f2 =
            √D[m] * α .* (bessely(n - 1, α * rD) - bessely(n + 1, α * rD)) / 2 +
            W[m] * bessely(n, α * rD)
    elseif dim == 3
        f1 =
            √D[m] * α .* (
                sphericalbesselj(n - 1, α * rD) - sphericalbesselj(n + 1, α * rD) -
                sphericalbesselj(n, α * rD) ./ (α * rD)
            ) / 2 + W[m] * sphericalbesselj(n, α * rD)
        f2 =
            √D[m] * α .* (
                sphericalbessely(n - 1, α * rD) - sphericalbessely(n + 1, α * rD) -
                sphericalbessely(n, α * rD) ./ (α * rD)
            ) / 2 + W[m] * sphericalbessely(n, α * rD)
    end

    bc1 .* f1 + bc2 .* f2
end

function compute_bc(params, α, n)
    (; r, D, W, dim) = params
    m = length(D)

    bc = zeros(m, 2)
    bc[:, 1] .= 1

    if α == 0
        return bc
    end

    ρ = r[1:(m-1)] ./ sqrt.(D[1:(m-1)])
    ρt = r[1:(m-1)] ./ sqrt.(D[2:m])

    if dim == 1
        bc[1, 1] = sqrt(D[1]) * α
    end

    for i = 1:(m-1)
        J, Y, dJ, dY = compute_JY(α * ρ[i], n, dim)
        J1, Y1, dJ1, dY1 = compute_JY(α * ρt[i], n, dim)

        A11 = -√(D[i] / D[i+1]) * dJ * Y1 + J * dY1 + √D[i] / W[i] * α * dJ * dY1
        A12 = -√(D[i] / D[i+1]) * dY * Y1 + Y * dY1 + √D[i] / W[i] * α * dY * dY1
        A21 = √(D[i] / D[i+1]) * dJ * J1 - J * dJ1 - √D[i] / W[i] * α * dJ * dJ1
        A22 = √(D[i] / D[i+1]) * dY * J1 - Y * dJ1 - √D[i] / W[i] * α * dY * dJ1

        z = α * ρt[i]
        if dim == 1
            Q = 1
        elseif dim == 2
            Q = z * π / 2
        elseif dim == 3
            Q = z^2
        end
        A = Q * [A11 A12; A21 A22]
        bc[i+1, :] = A * bc[i, :]
    end

    bc
end

function compute_β(params, α, n)
    I = compute_int_I(params, α, α, n)
    1 / √(2 * sum(I))
end

function compute_ν(n, dim)
    if dim == 1
        ν = 0
    elseif dim == 2
        ν = n^2
    elseif dim == 3
        ν = n * (n + 1)
    end
    ν
end

function compute_JY(z, n, dim)
    if dim == 1
        J = cos.(z)
        Y = sin.(z)
        dJ = -sin.(z)
        dY = cos.(z)
    elseif dim == 2
        J = besselj.(n, z)
        Y = bessely.(n, z)
        dJ = (besselj.(n - 1, z) - besselj.(n + 1, z)) / 2
        dY = (bessely.(n - 1, z) - bessely.(n + 1, z)) / 2
    elseif dim == 3
        J = sphericalbesselj.(n, z)
        Y = sphericalbessely.(n, z)

        dJ = [copy(z);]
        dY = [copy(z);]
        ind = z .> 0
        zz = z[ind]
        ind = [ind;]
        @. dJ[ind] =
            (
                sphericalbesselj(n - 1, zz) - sphericalbesselj(n + 1, zz) -
                sphericalbesselj(n, zz) / zz
            ) / 2
        @. dY[ind] =
            (
                sphericalbessely(n - 1, zz) - sphericalbessely(n + 1, zz) -
                sphericalbessely(n, zz) / zz
            ) / 2

        ind = [z .== 0;]
        dJ[ind] .= (n == 1) / 3
        dY[ind] .= Inf

        if size(z) == ()
            J, Y, dJ, dY = J[1], Y[1], dJ[1], dY[1]
        end
    end

    J, Y, dJ, dY
end

function compute_v(params, α, n, x, bc = compute_bc(params, α, n))
    (; r, D, dim) = params

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

function compute_int_I(params, α₁, α₂, n)
    (; r, D, dim) = params
    m = length(D)
    R = r[m]

    dr = r[m] * 1e-12
    r₂ = r[1:m] .- dr
    r₁ = [0; r[1:(m-1)]] .+ dr

    # Check if α₁ is not equal to α₂
    # if abs(α₁ - α₂) > 1e-8
    if !(abs(α₁ - α₂) ≈ 0)
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

function compute_int_K(params, α₁, n₁, α₂, n₂)
    (; r, D, dim) = params

    m = length(D)
    R = r[m]
    L = R

    dr = r[m] * 1e-12
    r₁ = [0; r[1:(m-1)]] .+ dr
    r₂ = r[1:m] .- dr

    K = 0
    if (dim == 1 && abs(α₁ - α₂) > 1e-6) || (dim > 1 && abs(n₁ - n₂) == 1)
        λ₁ = α₁^2
        λ₂ = α₂^2
        ν₁ = compute_ν(n₁, dim)
        ν₂ = compute_ν(n₂, dim)

        # Computation at point r₁
        V₁, dV₁ = compute_v(params, α₁, n₁, r₁)
        V₂, dV₂ = compute_v(params, α₂, n₂, r₁)
        I₁a = @.(
            ((λ₁ + λ₂) * D * r₁^(dim - 1) - (ν₁ + ν₂ - (dim - 1)) * D^2 * r₁^(dim - 3)) *
            V₁ *
            V₂
        )
        I₁b = @.(
            ((λ₂ - λ₁) * D * r₁^dim - (ν₁ - ν₂ - (dim - 1)) * D^2 * r₁^(dim - 2)) *
            dV₁ *
            V₂
        )
        I₁c = @.(
            ((λ₁ - λ₂) * D * r₁^dim - (ν₂ - ν₁ - (dim - 1)) * D^2 * r₁^(dim - 2)) *
            V₁ *
            dV₂
        )
        I₁d = @. 2 * D^2 * r₁^(dim - 1) * dV₁ * dV₂
        I₁ = I₁a + I₁b + I₁c + I₁d

        # Computation at point r₂
        V₁, dV₁ = compute_v(params, α₁, n₁, r₂)
        V₂, dV₂ = compute_v(params, α₂, n₂, r₂)
        I₂a = @.(
            ((λ₁ + λ₂) * D * r₂^(dim - 1) - (ν₁ + ν₂ - (dim - 1)) * D^2 * r₂^(dim - 3)) *
            V₁ *
            V₂
        )
        I₂b = @.(
            ((λ₂ - λ₁) * D * r₂^dim - (ν₁ - ν₂ - (dim - 1)) * D^2 * r₂^(dim - 2)) *
            dV₁ *
            V₂
        )
        I₂c = @.(
            ((λ₁ - λ₂) * D * r₂^dim - (ν₂ - ν₁ - (dim - 1)) * D^2 * r₂^(dim - 2)) *
            V₁ *
            dV₂
        )
        I₂d = @. 2 * D^2 * r₂^(dim - 1) * dV₁ * dV₂
        I₂ = I₂a + I₂b + I₂c + I₂d

        K = sum(I₂ - I₁) / (λ₁ - λ₂)^2 / R^(dim + 1)
    elseif dim == 1 && abs(α₁ - α₂) < 1e-6
        if α₁ > 1e-6
            λ = α₁^2
            V, dV = compute_v(params, α₁, n₁, r₁)
            I₁ = @. D * r₁^2 * dV^2 + (λ * r₁^2 + D) * V^2 - 2 * D * r₁ * dV * V

            V, dV = compute_v(params, α₁, n₁, r₂)
            I₂ = @. D * r₂^2 * dV^2 + (λ * r₂^2 + D) * V^2 - 2 * D * r₂ * dV * V

            K = sum(I₂ - I₁) / (4 * λ * L * R)
        else
            K = R / 2L
        end
    else
        K = 0
    end
    K
end
