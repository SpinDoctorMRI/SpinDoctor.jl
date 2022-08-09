"""
    solve(problem::AnalyticalLaplace)

Solve analytical Laplace eigenvalue problem for symmetrical geometries.
"""
function solve(laplace::AnalyticalLaplace)
    (; D, T₂, dim, eiglim, eigstep) = laplace

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

    (; Λ, β, U, B, Bri, Br, α, n)
end
