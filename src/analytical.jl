# Analytical matrix formalism for multilayered plates, cylinders and spheres.
#
# Based on the following articles and corresponding MATLAB code:
#   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
#       Rev. Mod. Phys. 79, 1077 (2007)
#   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted
#       diffusion in multilayered structures,
#       J. Magn. Reson. 205, 181-195 (2010).

"""
    solve_analytical_laplace(; rho, r, D, kappa, T2, gamma, dim, eiglim, eigstep)

Solve the analytical radial Laplace eigenvalue problem of multilayered plates (`dim = 1`),
cylinders (`dim = 2`) or spheres (`dim = 3`), with the coefficients of
[`analytical_coefficients`](@ref), for eigenvalues up to `eiglim`, searched with steps
`eigstep`. Return the coefficients together with the eigendecomposition.
"""
function solve_analytical_laplace(; rho, r, D, kappa, T2, gamma, dim, eiglim, eigstep)
    laplace = (; rho, r, D, kappa, T2, gamma, dim, eiglim, eigstep)

    L = length(D)

    # Compute radial and angular eigenvalues (lambda and nu), with alpha^2 = lambda and n^2 = nu for
    # cylinders and n(n+1) = nu for spheres
    dalpha = sqrt(eigstep) * 1e3
    alpha_min = dalpha / 100
    alpha_max = sqrt(eiglim) * 1e3
    alpha, n = find_alpha(laplace, alpha_min, alpha_max, dalpha)

    N = length(alpha)
    lambda = alpha .^ 2
    Lambda = diagm(lambda)

    beta = zeros(N)
    J = zeros(N)
    for m = 1:N
        beta[m] = compute_beta(laplace, alpha[m], n[m])
        if n[m] == 0
            J[m] = compute_int_J(laplace, alpha[m])[1]
        end
    end

    if dim == 1
        coeff = sqrt(2)
    elseif dim == 2
        coeff = 2
    elseif dim == 3
        coeff = sqrt(6)
    end
    U = coeff * J .* beta

    B = zeros(N, N)
    for m1 = 1:N
        for m2 = m1:N
            if dim == 1
                K = compute_int_K(laplace, alpha[m1], n[m1], alpha[m2], n[m2])
                B[m1, m2] = 2 * beta[m1] * beta[m2] * K
                B[m2, m1] = B[m1, m2]
            elseif abs(n[m1] - n[m2]) == 1
                K = compute_int_K(laplace, alpha[m1], n[m1], alpha[m2], n[m2])
                if dim == 2
                    epsilon = sqrt(1 + (n[m1] == 0) + (n[m2] == 0))
                elseif dim == 3
                    epsilon = (n[m1] + n[m2] + 1) / sqrt((2 * n[m1] + 1) * (2 * n[m2] + 1))
                end
                B[m1, m2] = epsilon * beta[m1] * beta[m2] * K
                B[m2, m1] = B[m1, m2]
            end
        end
    end

    Bri = zeros(N, N, L)
    for m1 = 1:N
        for m2 = m1:N
            if dim == 1
                int_I = compute_int_I(laplace, alpha[m1], alpha[m2], 0)
                Bri[m1, m2, :] = 2 * beta[m1] * beta[m2] * int_I
                Bri[m2, m1, :] = Bri[m1, m2, :]
            elseif n[m1] == n[m2]
                int_I = compute_int_I(laplace, alpha[m1], alpha[m2], n[m1])
                Bri[m1, m2, :] = 2 * beta[m1] * beta[m2] * int_I
                Bri[m2, m1, :] = Bri[m1, m2, :]
            end
        end
    end

    Br = zeros(N, N)
    for i = 1:L
        Br .+= Bri[:, :, i] / T2[i]
    end

    merge(laplace, (; Lambda, beta, U, B, Bri, Br, alpha, n))
end

"""
    solve_analytical_mf(laplace, volumes, gradient)

Compute the signal in multilayered plates, cylinders or spheres using an analytical matrix
formalism solution, from the eigendecomposition `laplace` of
[`solve_analytical_laplace`](@ref) and the compartment `volumes`. Only implemented for
[`PGSE`](@ref) sequences.

See [Grebenkov2007](@cite) and [Grebenkov2010a](@cite).
"""
function solve_analytical_mf(laplace, volumes, gradient::ScalarGradient)
    (; r, rho, gamma, Lambda, Br, B, U) = laplace

    gradient.profile isa PGSE || error("Only implemented for PGSE")

    w = gamma * gradient.amplitude * 1e12 * r[end]
    delta = gradient.profile.delta * 1e-6
    Delta = gradient.profile.Delta * 1e-6

    U0 = complex(U)
    U = complex(U)

    expmv!(-delta, Lambda + Br + im * w * B, U)
    expmv!(-(Delta - delta), Lambda + Br, U)
    expmv!(-delta, Lambda + Br - im * w * B, U)

    U0'U * volumes'rho
end

function find_alpha(params, alpha_min, alpha_max, dalpha)
    (; dim) = params

    if dim == 1
        alpha = find_zeros(alpha -> alpha_func(params, alpha, 0), alpha_min, alpha_max)
        n = zeros(Int, size(alpha))
    else
        alpha = Float64[]
        n = Int[]
        alpha0 = [0.0]
        j = 0
        while !isempty(alpha0) && alpha0[1] < alpha_max
            @debug "Computing zeros of F_n(alpha) for n = $j"
            alpha0 = find_zeros(alpha -> alpha_func(params, alpha, j), alpha_min, alpha_max)
            append!(alpha, alpha0)
            append!(n, fill(j, size(alpha0)))
            j = j + 1
        end

        inds = sortperm(alpha)
        alpha = alpha[inds]
        n = n[inds]
    end

    # Number of layers
    L = length(params.D)

    inds_keep = dalpha / 100 .< alpha .< alpha_max
    alpha = alpha[inds_keep]
    n = n[inds_keep]

    if params.kappa[L] < 1e-12 && (length(params.kappa) == L || params.kappa[end] < 1e-12)
        # Add ground eigenmode explicitly
        pushfirst!(alpha, 0)
        pushfirst!(n, 0)
    end

    alpha, n
end

function alpha_func(params, alpha, n)
    (; r, D, kappa, dim) = params

    m = length(D)
    bc1 = zeros(size(alpha))
    bc2 = zeros(size(alpha))

    for i = 1:length(alpha)
        bc = compute_bc(params, alpha[i], n)
        bc1[i] = bc[m, 1]
        bc2[i] = bc[m, 2]
    end

    rD = r[m] / sqrt(D[m])

    if dim == 1
        f1 = -sqrt(D[m]) * alpha .* sin(alpha * rD) + kappa[m] * cos(alpha * rD)
        f2 = sqrt(D[m]) * alpha .* cos(alpha * rD) + kappa[m] * sin(alpha * rD)
    elseif dim == 2
        f1 =
            sqrt(D[m]) * alpha .*
            (besselj(n - 1, alpha * rD) - besselj(n + 1, alpha * rD)) / 2 +
            kappa[m] * besselj(n, alpha * rD)
        f2 =
            sqrt(D[m]) * alpha .*
            (bessely(n - 1, alpha * rD) - bessely(n + 1, alpha * rD)) / 2 +
            kappa[m] * bessely(n, alpha * rD)
    elseif dim == 3
        f1 =
            sqrt(D[m]) * alpha .* (
                sphericalbesselj(n - 1, alpha * rD) - sphericalbesselj(n + 1, alpha * rD) -
                sphericalbesselj(n, alpha * rD) ./ (alpha * rD)
            ) / 2 + kappa[m] * sphericalbesselj(n, alpha * rD)
        f2 =
            sqrt(D[m]) * alpha .* (
                sphericalbessely(n - 1, alpha * rD) - sphericalbessely(n + 1, alpha * rD) -
                sphericalbessely(n, alpha * rD) ./ (alpha * rD)
            ) / 2 + kappa[m] * sphericalbessely(n, alpha * rD)
    end

    bc1 .* f1 + bc2 .* f2
end

function compute_bc(params, alpha, n)
    (; r, D, kappa, dim) = params
    m = length(D)

    bc = zeros(m, 2)
    bc[:, 1] .= 1

    if alpha == 0
        return bc
    end

    rho = r[1:(m-1)] ./ sqrt.(D[1:(m-1)])
    rho_t = r[1:(m-1)] ./ sqrt.(D[2:m])

    if dim == 1
        bc[1, 1] = sqrt(D[1]) * alpha
    end

    for i = 1:(m-1)
        J, Y, dJ, dY = compute_JY(alpha * rho[i], n, dim)
        J1, Y1, dJ1, dY1 = compute_JY(alpha * rho_t[i], n, dim)

        A11 =
            -sqrt(D[i] / D[i+1]) * dJ * Y1 +
            J * dY1 +
            sqrt(D[i]) / kappa[i] * alpha * dJ * dY1
        A12 =
            -sqrt(D[i] / D[i+1]) * dY * Y1 +
            Y * dY1 +
            sqrt(D[i]) / kappa[i] * alpha * dY * dY1
        A21 =
            sqrt(D[i] / D[i+1]) * dJ * J1 - J * dJ1 -
            sqrt(D[i]) / kappa[i] * alpha * dJ * dJ1
        A22 =
            sqrt(D[i] / D[i+1]) * dY * J1 - Y * dJ1 -
            sqrt(D[i]) / kappa[i] * alpha * dY * dJ1

        z = alpha * rho_t[i]
        if dim == 1
            Q = 1
        elseif dim == 2
            Q = z * pi / 2
        elseif dim == 3
            Q = z^2
        end
        A = Q * [A11 A12; A21 A22]
        bc[i+1, :] = A * bc[i, :]
    end

    bc
end

function compute_beta(params, alpha, n)
    I = compute_int_I(params, alpha, alpha, n)
    1 / sqrt(2 * sum(I))
end

function compute_nu(n, dim)
    if dim == 1
        nu = 0
    elseif dim == 2
        nu = n^2
    elseif dim == 3
        nu = n * (n + 1)
    end
    nu
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

function compute_v(params, alpha, n, x, bc = compute_bc(params, alpha, n))
    (; r, D, dim) = params

    x = [x;]

    m = length(D)
    v = ones(size(x))
    dv = zeros(size(x))
    if alpha > 0
        r1 = [0; r[1:m]]
        r1[1] -= 1e-8
        r1[m+1] += 1e-8
        for i = 1:m
            inds = r1[i] .< x .<= r1[i+1]
            J, Y, dJ, dY = compute_JY(alpha * x[inds] / sqrt(D[i]), n, dim)
            if abs(bc[i, 2]) > 1e-12
                v[inds] = bc[i, 1] * J + bc[i, 2] * Y
                dv[inds] = alpha / sqrt(D[i]) * (bc[i, 1] * dJ + bc[i, 2] * dY)
            else
                v[inds] = bc[i, 1] * J
                dv[inds] = alpha / sqrt(D[i]) * bc[i, 1] * dJ
            end
        end
    end

    v, dv
end

function compute_int_I(params, alpha1, alpha2, n)
    (; r, D, dim) = params
    m = length(D)
    R = r[m]

    dr = r[m] * 1e-12
    r2 = r[1:m] .- dr
    r1 = [0; r[1:(m-1)]] .+ dr

    # Check if alpha1 is not equal to alpha2
    if !iszero(alpha1 - alpha2)
        lambda1 = alpha1^2
        lambda2 = alpha2^2

        # Computation at point r1
        V1, dV1 = compute_v(params, alpha1, n, r1)
        V2, dV2 = compute_v(params, alpha2, n, r1)
        I1 = D .* r1 .^ (dim - 1) .* (dV1 .* V2 - dV2 .* V1)

        # Computation at point r2
        V1, dV1 = compute_v(params, alpha1, n, r2)
        V2, dV2 = compute_v(params, alpha2, n, r2)
        I2 = D .* r2 .^ (dim - 1) .* (dV1 .* V2 - dV2 .* V1)

        I = (I2 - I1) / (lambda1 - lambda2) / R^dim
    else
        alpha = alpha1
        lambda = alpha^2
        nu = compute_nu(n, dim)

        if alpha > 1e-6
            V, dV = compute_v(params, alpha, n, r1)
            I1 = (
                D .* r1 .^ dim .* dV .^ 2 +
                (lambda * r1 .^ dim - D .* r1 .^ (dim - 2) * nu) .* V .^ 2 +
                (dim - 2) * D .* r1 .^ (dim - 1) .* dV .* V
            )

            V, dV = compute_v(params, alpha, n, r2)
            I2 = (
                D .* r2 .^ dim .* dV .^ 2 +
                (lambda * r2 .^ dim - D .* r2 .^ (dim - 2) * nu) .* V .^ 2 +
                (dim - 2) * D .* r2 .^ (dim - 1) .* dV .* V
            )

            I = (I2 - I1) / (2 * lambda * R^dim)
        else
            I = (r2 .^ dim - r1 .^ dim) / R^dim / dim
        end
    end

    I
end

function compute_int_J(params, alpha)
    (; r, D, kappa, dim) = params
    m = length(D)
    R = r[m]

    lambda = alpha^2

    if alpha > 1e-6
        V, _ = compute_v(params, alpha, 0, R)
        J = kappa[m] * V / lambda / R
    else
        J = 1 / dim
    end

    J
end

function compute_int_K(params, alpha1, n1, alpha2, n2)
    (; r, D, dim) = params

    m = length(D)
    R = r[m]
    L = R

    dr = r[m] * 1e-12
    r1 = [0; r[1:(m-1)]] .+ dr
    r2 = r[1:m] .- dr

    K = 0
    if (dim == 1 && abs(alpha1 - alpha2) > 1e-6) || (dim > 1 && abs(n1 - n2) == 1)
        lambda1 = alpha1^2
        lambda2 = alpha2^2
        nu1 = compute_nu(n1, dim)
        nu2 = compute_nu(n2, dim)

        # Computation at point r1
        V1, dV1 = compute_v(params, alpha1, n1, r1)
        V2, dV2 = compute_v(params, alpha2, n2, r1)
        I1a = @.(
            (
                (lambda1 + lambda2) * D * r1^(dim - 1) -
                (nu1 + nu2 - (dim - 1)) * D^2 * r1^(dim - 3)
            ) *
            V1 *
            V2
        )
        I1b = @.(
            (
                (lambda2 - lambda1) * D * r1^dim -
                (nu1 - nu2 - (dim - 1)) * D^2 * r1^(dim - 2)
            ) *
            dV1 *
            V2
        )
        I1c = @.(
            (
                (lambda1 - lambda2) * D * r1^dim -
                (nu2 - nu1 - (dim - 1)) * D^2 * r1^(dim - 2)
            ) *
            V1 *
            dV2
        )
        I1d = @. 2 * D^2 * r1^(dim - 1) * dV1 * dV2
        I1 = I1a + I1b + I1c + I1d

        # Computation at point r2
        V1, dV1 = compute_v(params, alpha1, n1, r2)
        V2, dV2 = compute_v(params, alpha2, n2, r2)
        I2a = @.(
            (
                (lambda1 + lambda2) * D * r2^(dim - 1) -
                (nu1 + nu2 - (dim - 1)) * D^2 * r2^(dim - 3)
            ) *
            V1 *
            V2
        )
        I2b = @.(
            (
                (lambda2 - lambda1) * D * r2^dim -
                (nu1 - nu2 - (dim - 1)) * D^2 * r2^(dim - 2)
            ) *
            dV1 *
            V2
        )
        I2c = @.(
            (
                (lambda1 - lambda2) * D * r2^dim -
                (nu2 - nu1 - (dim - 1)) * D^2 * r2^(dim - 2)
            ) *
            V1 *
            dV2
        )
        I2d = @. 2 * D^2 * r2^(dim - 1) * dV1 * dV2
        I2 = I2a + I2b + I2c + I2d

        K = sum(I2 - I1) / (lambda1 - lambda2)^2 / R^(dim + 1)
    elseif dim == 1 && abs(alpha1 - alpha2) < 1e-6
        if alpha1 > 1e-6
            lambda = alpha1^2
            V, dV = compute_v(params, alpha1, n1, r1)
            I1 = @. D * r1^2 * dV^2 + (lambda * r1^2 + D) * V^2 - 2 * D * r1 * dV * V

            V, dV = compute_v(params, alpha1, n1, r2)
            I2 = @. D * r2^2 * dV^2 + (lambda * r2^2 + D) * V^2 - 2 * D * r2 * dV * V

            K = sum(I2 - I1) / (4 * lambda * L * R)
        else
            K = R / 2L
        end
    else
        K = 0
    end
    K
end
