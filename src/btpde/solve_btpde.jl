"""
    solve_btpde(model, experiment)

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
"""
function solve_btpde(model::Model, experiment::Experiment)

    # Measure function evalutation time
    starttime = Base.time()

    # Extract input parameters
    @unpack mesh, D, T₂, ρ = model
    @unpack directions, sequences, values, values_type = experiment.gradient
    @unpack odesolver, reltol, abstol, nsave = experiment.btpde

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
    Mx_cmpts = [[] for _ = 1:3]
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
            push!(Mx_cmpts[dim], assemble_mass_matrix(elements', volumes, points[dim, :]))
        end
    end

    # Assemble global finite element matrices
    M = blockdiag(M_cmpts...)
    S = blockdiag(S_cmpts...)
    R = blockdiag((M_cmpts ./ T₂)...)
    Mx = [blockdiag(Mx_cmpts[dim]...) for dim = 1:3]
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

    # Time dependent ODE function
    function Mdξ!(dξ, ξ, p, t)
        # @show t
        @unpack J, S, Q, A, R, q, f = p
        @. J = -(S + Q + R + im * f(t) * q * A)
        mul!(dξ, J, ξ)
        nothing
    end

    # Time independent ODE function, given Jacobian
    function Mdξ_constant!(dξ, ξ, p, t)
        # @show t
        J = p.J
        mul!(dξ, J, ξ)
        nothing
    end

    # Time dependent Jacobian of ODE function with respect to the state `ξ`
    function Jac!(J, ξ, p, t)
        @unpack S, Q, A, R, q, f = p
        @. J = -(S + Q + R + im * f(t) * q * A)
        nothing
    end

    # Jacobian sparsity pattern
    jac_prototype = -(S + Q + im * sum(Mx))

    # Iterate over gradient amplitudes, time profiles and directions
    for idir = 1:ndirection, iseq = 1:nsequence, iamp = 1:namplitude

        # Gradient amplitude
        q = qvalues[iamp, iseq]
        b = bvalues[iamp, iseq]

        # Time profile
        f = sequences[iseq]
        interval = (0, echotime(f))
        if nsave == 1
            saveat = [interval[2]]
            # saveat = [interval[1], interval[2]]
        else
            saveat = LinRange(interval..., nsave)
        end

        # Gradient direction
        dir = directions[:, idir]

        # Display state of iterations
        @printf "Solving BTPDE with size %d\n" (sum(npoint_cmpts))
        @printf "  Direction %d of %d: g = [%.2f, %.2f, %.2f]\n" idir ndirection dir...
        @printf "  Sequence  %d of %d: f = %s\n" iseq nsequence f
        @printf "  Amplitude %d of %d: q = %g, b = %g\n" iamp namplitude q b

        # Gradient direction dependent finite element matrix
        A = dir' * Mx

        # ODE problem
        J = -(S + Q + im * q * A)
        p = (; J, S, Q, A, R, q, f)

        # Gather ODE function
        if all(constant_intervals(p.f)) # is_constant(p.f, t)
            func = Mdξ_constant!
            Jac!(p.J, ρ, p, 0)
            # jac = (J, _, p, t) -> (@show t; J .= p.J; nothing)
            jac = (J, _, p, t) -> (J .= p.J; nothing)
        else
            func = Mdξ!
            jac = Jac!
        end
        odefunction =
            ODEFunction(func, jac = jac, jac_prototype = jac_prototype, mass_matrix = M)
        odeproblem = ODEProblem(odefunction, ρ, interval, p, progress = false)#, dtmax = 50)

        tstops = intervals(f)[2:end-1]

        # Callback for updating problem
        function affect!(integrator)
            @unpack u, p, t = integrator
            # print("    t = $t: constant time profile, f = $(p.f(t))")
            Jac!(p.J, u, p, t)
        end
        callback = PresetTimeCallback(tstops, affect!)

        # Measure solving time
        itertime = Base.time()

        # Solve ODE problem
        sol = solve(
            odeproblem,
            odesolver,
            saveat = saveat,
            reltol = reltol,
            abstol = abstol,
            callback = callback,
        )

        itertimes[iamp, iseq, idir] = Base.time() - itertime

        # Extract solution
        time[iamp, iseq, idir] = sol.t
        ξ = hcat(sol.u...)

        if nsave == 1
            time[iamp, iseq, idir] =  time[iamp, iseq, idir][[end]]
            ξ = ξ[:, [end]]
        end

        # Split solution into compartments
        for icmpt = 1:ncompartment
            inds = inds_cmpts[icmpt]+1:inds_cmpts[icmpt+1]

            # Store magnetization in compartment
            magnetization[icmpt, iamp, iseq, idir] = ξ[inds, :]

            # Integrate final magnetization over compartment
            signal[icmpt, iamp, iseq, idir] =
                sum(M_cmpts[icmpt] * ξ[inds, end], dims = 1)[1]

        end
    end

    signal_allcmpts = sum(signal, dims = 1)[1, :, :, :]

    totaltime = Base.time() - starttime

    (; magnetization, signal, signal_allcmpts, time, itertimes, totaltime)
end
