"""
    solve_btpde_midpoint(model, setup)

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
This function uses a manual time stepping scheme (theta-rule), that requires a degree of
implicitness `θ` and a time step `dt`.
    `θ = 0.5`: Crank-Nicolson (second order)
    `θ = 1.0`: Implicit Euler (first order)
"""
function solve_btpde_midpoint(model::Model, experiment::Experiment)

    # Measure function evalutation time
    starttime = Base.time()

    # Extract input parameters
    @unpack mesh, D, T₂, ρ = model
    @unpack directions, sequences, values, values_type = experiment.gradient
    θ = experiment.btpde_midpoint.θ
    dt = experiment.btpde_midpoint.timestep

    # Deduce sizes
    ncompartment = length(ρ)
    npoint_cmpts = size.(mesh.points, 2)
    inds_cmpts = cumsum([0; npoint_cmpts])
    ndirection = size(directions, 2)
    nsequence = length(sequences)
    namplitude = length(values)

    # Assemble finite element matrices compartment-wise
    M_cmpts = []
    S_cmpts = []
    Ax_cmpts = [[] for _ = 1:3]
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        volumes, _ = get_mesh_volumes(points, elements)

        # Assemble mass, stiffness and flux matrices
        push!(M_cmpts, assemble_mass_matrix(elements', volumes))
        push!(S_cmpts, assemble_stiffness_matrix(elements', points', D[icmpt]))

        # Assemble first order product moment matrices
        for dim = 1:3
            push!(Ax_cmpts[dim], assemble_mass_matrix(elements', volumes, points[dim, :]))
        end
    end

    # Assemble global finite element matrices
    M = blockdiag(M_cmpts...)
    S = blockdiag(S_cmpts...)
    R = blockdiag((M_cmpts ./ T₂)...)
    Ax = [blockdiag(Ax_cmpts[dim]...) for dim = 1:3]
    Q_blocks = assemble_flux_matrices(mesh.points, mesh.facets)
    Q = couple_flux_matrix(model, Q_blocks, false)

    # Create initial conditions (enforce complex values)
    ρ = vcat(fill.(complex.(ρ), npoint_cmpts)...)

    # Allocate output arrays
    signal = zeros(ComplexF64, ncompartment, namplitude, nsequence, ndirection)
    signal_allcmpts = zeros(ComplexF64, namplitude, nsequence, ndirection)
    magnetization =
        Array{Matrix{ComplexF64},4}(undef, ncompartment, namplitude, nsequence, ndirection)
    time = Array{Vector{Float64},3}(undef, namplitude, nsequence, ndirection)
    itertimes = zeros(namplitude, nsequence, ndirection)

    # Q-values and b-values
    if values_type == "q"
        qvalues = repeat(values, 1, nsequence)
        bvalues = values .^ 2 .* bvalue_no_q.(sequences)'
    else
        bvalues = repeat(values, 1, nsequence)
        qvalues = .√(values ./ bvalue_no_q.(sequences)')
    end

    Jac(q, fₜ, g) = -(S + Q + R + im * fₜ * q * (g' * Ax))

    # Iterate over gradient amplitudes, time profiles and directions
    for idir = 1:ndirection, iseq = 1:nsequence, iamp = 1:namplitude

        # Measure solving time
        itertime = Base.time()

        # Gradient direction
        dir = directions[:, idir]

        # Time profile
        f = sequences[iseq]
        ivals = intervals(f)

        # Gradient amplitude
        q = qvalues[iamp, iseq]
        b = bvalues[iamp, iseq]

        # Display state of iterations
        @printf "Solving BTPDE with size %d\n" (sum(npoint_cmpts))
        @printf "  Direction %d of %d: g = [%.2f, %.2f, %.2f]\n" idir ndirection dir...
        @printf "  Sequence  %d of %d: f = %s\n" iseq nsequence f
        @printf "  Amplitude %d of %d: q = %g, b = %g\n" iamp namplitude q b

        # Crank-Nicolson time stepping
        t = 0
        y = copy(ρ)
        Ey = copy(y)
        for i = 1:length(ivals)-1
            @printf "    Solving for interval [%g, %g]\n" ivals[i] ivals[i+1]
            J = Jac(q, f(ivals[i]), dir)
            F = lu(M - dt * θ * J)
            # F = factorize(M - dt * θ * J)
            E = M + dt * (1 - θ) * J
            while t + dt < ivals[i+1]
                # @show t
                mul!(Ey, E, y)
                ldiv!(y, F, Ey)
                t += dt
            end
            dt_last = ivals[i+1] - t
            mul!(Ey, M + dt_last * (1 - θ) * J, y)
            ldiv!(y, lu(M - dt_last * θ * J), Ey)
            t += dt_last
        end

        # Extract solution
        time[iamp, iseq, idir] = [t]
        mag = [y;]

        # Split solution into compartments
        for icmpt = 1:ncompartment
            inds = inds_cmpts[icmpt]+1:inds_cmpts[icmpt+1]

            # Store magnetization in compartment
            magnetization[icmpt, iamp, iseq, idir] = mag[inds, :]

            # Integrate final magnetization over compartment
            signal[icmpt, iamp, iseq, idir] =
                sum(M_cmpts[icmpt] * mag[inds, end], dims = 1)[1]
        end

        itertimes[iamp, iseq, idir] = Base.time() - itertime
    end

    signal_allcmpts = sum(signal, dims = 1)[1, :, :, :]

    totaltime = Base.time() - starttime

    (; magnetization, signal, signal_allcmpts, time, itertimes, totaltime)
end
