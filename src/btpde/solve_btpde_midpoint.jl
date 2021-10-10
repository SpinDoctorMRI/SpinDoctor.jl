"""
    solve_btpde_midpoint(model, matrices, experiment)

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
This function uses a manual time stepping scheme (theta-rule), that requires a degree of
implicitness `θ` and a time step `Δt`.
    `θ = 0.5`: Crank-Nicolson (second order)
    `θ = 1.0`: Implicit Euler (first order)
"""
function solve_btpde_midpoint(model::Model, matrices, experiment::Experiment)

    # Measure function evalutation time
    starttime = Base.time()

    # Extract input parameters
    @unpack mesh, D, T₂, ρ = model
    @unpack M, S, R, Mx, Q, M_cmpts = matrices
    @unpack directions, sequences, values, values_type = experiment.gradient
    @unpack θ, timestep = experiment.btpde_midpoint

    # Deduce sizes
    ncompartment = length(ρ)
    npoint_cmpts = size.(mesh.points, 2)
    inds_cmpts = cumsum([0; npoint_cmpts])
    ndirection = size(directions, 2)
    nsequence = length(sequences)
    namplitude = length(values)

    qvalues, bvalues = get_values(experiment.gradient)

    # Create initial conditions (enforce complex values)
    ρ = vcat(fill.(complex.(ρ), npoint_cmpts)...)

    # Allocate output arrays
    signal = zeros(ComplexF64, ncompartment, namplitude, nsequence, ndirection)
    signal_allcmpts = zeros(ComplexF64, namplitude, nsequence, ndirection)
    magnetization =
        Array{Matrix{ComplexF64},4}(undef, ncompartment, namplitude, nsequence, ndirection)
    time = Array{Vector{Float64},3}(undef, namplitude, nsequence, ndirection)
    itertimes = zeros(namplitude, nsequence, ndirection)

    Jac(q, fₜ, g) = -(S + Q + R + im * fₜ * q * (g' * Mx))

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
        t = 0.0
        y = copy(ρ)
        Ey = copy(y)
        for i = 1:length(ivals)-1
            @printf "    Solving for interval [%g, %g]\n" ivals[i] ivals[i+1]

            # Adjust time step to obtain divide interval uniformly
            ival_length = ivals[i+1] - ivals[i]
            nt = round(Int, ival_length / timestep)
            Δt = ival_length / nt

            # Build matrices for interval
            J = Jac(q, f(ivals[i]), dir)
            F = lu(M .- Δt .* θ .* J)
            # F = factorize(M .- Δt .* θ .* J)
            E = @. complex(M) + Δt * (1 - θ) * J

            # Advance one interval
            # for it = 1:nt
            while t + Δt ≤ ivals[i+1]
                # @show t
                mul!(Ey, E, y)
                ldiv!(y, F, Ey)
                t += Δt
            end
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
