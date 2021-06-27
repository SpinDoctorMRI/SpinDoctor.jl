"""
    solve_analytical(setup, experiment[, volumes])

Compute the signal in a multilayered cylinder or sphere using an analytical
matrix formalism solution.

This function is based on the following articles and corresponding code:
    [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
        Rev. Mod.Phys. 79, 1077 (2007)
    [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted diffusion in
        multilayered structures, J. Magn. Reson. 205, 181-195 (2010).
"""
function solve_analytical(setup, experiment, volumes = nothing)

    # Sizes
    namplitude = length(experiment.gradient.values)
    nsequence = length(experiment.gradient.sequences)

    # Extract parameters
    sequences = experiment.gradient.sequences
    rmean = (setup.rmin + setup.rmax) / 2

    # Q-values and b-values
    if experiment.gradient.values_type == "q"
        qvalues = repeat(experiment.gradient.values, 1, nsequence)
        bvalues =
            experiment.gradient.values .^ 2 .* bvalue_no_q.(experiment.gradient.sequences)'
    else
        bvalues = repeat(experiment.gradient.values, 1, nsequence)
        qvalues =
            .√(experiment.gradient.values ./ bvalue_no_q.(experiment.gradient.sequences)')
    end

    if isa(setup, CylinderSetup)
        dim = 2
        get_vol = r -> π * r^2 * setup.height
    else
        dim = 3
        get_vol = r -> 4 * π / 3 * r^3
    end

    # Get OUT parameters
    ρ_out = setup.ρ_out
    r_out = (setup.rmin + setup.rmax) / 2
    D_out = tr(setup.D_out) / 3
    T_out = setup.T₂_out
    W_out = 0

    # Get IN parameters
    if setup.include_in
        ρ_in = setup.ρ_in
        r_in = setup.in_ratio * r_out
        D_in = tr(setup.D_in) / 3
        W_in = setup.κ_in_out
        T_in = setup.T₂_in
    else
        # Empty IN
        ρ_in = []
        r_in = []
        D_in = []
        W_in = []
        T_in = []
    end

    # Get ECS parameters
    if setup.ecs_shape != "no_ecs"
        # Include ECS
        ρ_ecs = setup.ρ_ecs
        r_ecs = r_out + setup.ecs_ratio * rmean
        D_ecs = tr(setup.D_ecs) / 3
        W_ecs = setup.κ_ecs
        T_ecs = setup.T₂_ecs

        # Add interface permeability between out and ecs
        W_out = setup.κ_out_ecs
    else
        # Empty ECS
        ρ_ecs = []
        r_ecs = []
        D_ecs = []
        W_ecs = []
        T_ecs = []

        if dim == 3
            # Add surface relaxivity for the outermost sphere ("out")
            W_out = setup.κ_out
        end
    end

    # Create parameters structure
    ρ = vcat(ρ_in, ρ_out, ρ_ecs)
    r = vcat(r_in, r_out, r_ecs) * 1e-6
    D = vcat(D_in, D_out, D_ecs) * 1e-6
    W = vcat(W_in, W_out, W_ecs) * 1
    T = vcat(T_in, T_out, T_ecs) * 1e-6
    params = (; ρ, r, D, W, T, dim)

    L = length(D)

    # Get spherical volumes
    volumes_exact = get_vol.(r * 1e6)

    # Subtract inner volumes (each volume is contained within the next)
    volumes_exact[2:end] = volumes_exact[2:end] - volumes_exact[1:end-1]


    # Volumes
    if isnothing(volumes)
        display("Using exact volumes")
        volumes = volumes_exact
    else
        display("Exact volumes:")
        display(volumes_exact)
        display("Finite element volumes:")
        display(volumes)
        display("Relative error:")
        display(abs.(volumes - volumes_exact) ./ volumes_exact)
        display("Using finite element volumes")
    end

    # Volume weighted quantities
    D_mean = D' * volumes / sum(volumes)
    S0 = volumes' * ρ

    # Upper Laplace eigenvalue limit
    eiglim = π^2 * D_mean / experiment.analytical.length_scale^2

    # Compute radial and angular eigenvalues (λ and ν), with
    # α^2 = λ and n^2 = ν for cylinders and n(n+1) = ν for spheres
    dα = √experiment.analytical.eigstep * 1e3
    α_min = dα / 100
    α_max = √eiglim * 1e6
    α, n = find_α(params, α_min, α_max, dα)

    N = length(α)
    λ = α .^ 2
    Λ = diagm(λ)

    β = zeros(N)
    J = zeros(N)
    for m = 1:N
        β[m] = compute_β(params, α[m], n[m])
        if n[m] == 0
            J[m] = compute_int_J(params, α[m])[1]
        end
    end

    if dim == 2
        coeff = 2
    else
        coeff = √6
    end
    U = coeff * J .* β

    B = zeros(N, N)
    for m1 = 1:N
        for m2 = m1:N
            if abs(n[m1] - n[m2]) == 1
                K = compute_int_K(params, α[m1], n[m1], α[m2], n[m2])
                if dim == 2
                    ε = √(1 + (n[m1] == 0) + (n[m2] == 0))
                else
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
            if n[m1] == n[m2]
                int_I = compute_int_I(params, α[m1], α[m2], n[m1])
                Bri[m1, m2, :] = 2 * β[m1] * β[m2] * int_I
                Bri[m2, m1, :] = Bri[m1, m2, :]
            end
        end
    end

    Br = zeros(N, N)
    for i = 1:L
        Br .+= Bri[:, :, i] / T[i]
    end

    # Iterate over experiments
    signal = zeros(namplitude, nsequence)
    for iseq = 1:nsequence, iamp = 1:namplitude
        seq = sequences[iseq]
        w = qvalues[iamp, iseq] * 1e12 * r[end]

        δ = seq.δ * 1e-6
        Δ = seq.Δ * 1e-6

        # Compute signal attenuation
        E = (
            U' *
            exp(-δ * (Λ + Br + im * w * B)) *
            exp(-(Δ - δ) * (Λ + Br)) *
            exp(-δ * (Λ + Br - im * w * B)) *
            U
        )
        signal[iamp, iseq] = E * S0
    end

    signal
end
