"""
    results = solve_btpde(mesh, setup)

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
"""
function solve_btpde(mesh, setup)

    # Measure time of function evaluation
    starttime = Base.time()

    # Extract input parameters
    @unpack boundary_markers, compartments, boundaries, σ, T₂, κ, ρ = setup.pde
    @unpack directions, sequences, values, values_type = setup.gradient
    @unpack odesolver, reltol, abstol, nsave = setup.btpde

    # Deduce sizes
    ncompartment, nboundary = size(boundary_markers)
    npoint_cmpts = size.(mesh.points, 2)
    inds_cmpts = cumsum([0; npoint_cmpts])
    ndirection = size(directions, 2)
    nsequence = length(sequences)
    namplitude = length(values)

    # Assemble finite element matrices compartment-wise
    fem_mat_cmpts = (M = [], S = [], Q = [], Mx = [[] for idim = 1:3])
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt]
        facets = mesh.facets[icmpt, :]
        elements = mesh.elements[icmpt]
        volumes, _ = get_mesh_volumes(points, elements)

        # Assemble mass, stiffness and flux matrices
        push!(fem_mat_cmpts.M, assemble_mass_matrix(elements', volumes))
        push!(fem_mat_cmpts.S, assemble_stiffness_matrix(elements', points', σ[icmpt]))
        push!(fem_mat_cmpts.Q, assemble_flux_matrix_cmpt(points, facets, κ))

        # Assemble first order product moment matrices
        for dim = 1:3
            push!(
                fem_mat_cmpts.Mx[dim],
                assemble_mass_matrix(elements', volumes, points[dim, :]),
            )
        end
    end

    # Assemble global finite element matrices
    M = blockdiag(fem_mat_cmpts.M...)
    S = blockdiag(fem_mat_cmpts.S...)
    Q = couple_flux_matrix(mesh, fem_mat_cmpts.Q)
    Mx = [blockdiag(fem_mat_cmpts.Mx[dim]...) for dim = 1:3]
    R = blockdiag((fem_mat_cmpts.M ./ T₂)...)

    # # Display sparsity
    # println("Mass matrix:")
    # display(spy(M))
    # println("Stiffness matrix:")
    # display(spy(S))
    # println("Flux matrix:")
    # display(spy(Q))

    # Create initial conditions (enforce complex values)
    ρ = vcat(fill.(complex(ρ), npoint_cmpts)...)

    # Allocate output arrays
    signal = zeros(ComplexF64, ncompartment, namplitude, nsequence, ndirection)
    signal_allcmpts = zeros(ComplexF64, namplitude, nsequence, ndirection)
    magnetization = fill(
        Array{ComplexF64,2}(undef, 0, 0),
        ncompartment,
        namplitude,
        nsequence,
        ndirection,
    )
    time = fill(Float64[], namplitude, nsequence, ndirection)
    itertimes = zeros(namplitude, nsequence, ndirection)

    # Q-values and b-values
    if values_type == 'q'
        qvalues = repeat(values, 1, nsequence)
        bvalues = values .^ 2 .* bvalue_no_q.(sequences)'
    else
        bvalues = repeat(values, 1, nsequence)
        qvalues = .√(values ./ bvalue_no_q.(sequences)')
    end

    # Time dependent ODE function
    function M∂u∂t!(du, u, p, t)
        J, S, Q, A, R, q, f = p
        @. J = -(S + Q + R + im * f(t) * q * A)
        mul!(du, J, u)
    end

    # Time independent ODE function, given jacobian `J` 
    function M∂u∂t_constant!(du, u, p, t)
        J, _ = p
        mul!(du, J, u)
    end

    # Time dependent Jacobian of ODE function with respect to the state `u`
    function Jac!(J, u, p, t)
        _, S, Q, A, R, q, f = p
        @. J = -(S + Q + R + im * f(t) * q * A)
    end

    # Jacobian sparsity pattern
    jac_prototype = -(S + Q + im * sum(Mx))

    # Iterate over gradient amplitudes, time profiles and directions
    for iamp = 1:namplitude, iseq = 1:nsequence, idir = 1:ndirection

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
        A = sum(dir[dim] * Mx[dim] for dim = 1:3)

        # ODE problem
        J = -(S + Q + im * q * A)
        p = (; J, S, Q, A, R, q, f)

        # Gather ODE function
        function get_odefunction(u, t, p)
            print("    t = $t: ")
            if all(constant_intervals(p.f)) # is_constant(p.f, t)
                println("constant time profile, f = $(p.f(t))")
                func = M∂u∂t_constant!
                Jac!(p.J, u, p, t)
                jac = (J, u, p, t) -> (J .= p.J)
            else
                println("time dependent time profile")
                func = M∂u∂t!
                jac = Jac!
            end
            odefunction =
                ODEFunction(func, jac = jac, jac_prototype = jac_prototype, mass_matrix = M)
        end
        odeproblem = ODEProblem(get_odefunction(ρ, 0, p), ρ, interval, p, progress = false)

        tstops = intervals(f)[2:end-1]

        # Callback for updating problem
        function affect!(integrator)
            integrator.f = get_odefunction(integrator.u, integrator.t, integrator.p)
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
        mag = hcat(sol.u...)

        if nsave == 1
            mag = mag[:, [end]]
        end

        # Split solution into compartments
        for icmpt = 1:ncompartment
            inds = inds_cmpts[icmpt]+1:inds_cmpts[icmpt+1]

            # Store magnetization in compartment
            magnetization[icmpt, iamp, iseq, idir] = mag[inds, :]

            # Integrate final magnetization over compartment
            signal[icmpt, iamp, iseq, idir] =
                sum(fem_mat_cmpts.M[icmpt] * mag[inds, end], dims = 1)[1]

        end
    end

    signal_allcmpts = sum(signal, dims = 1)[1, :, :, :]

    totaltime = Base.time() - starttime

    (; magnetization, signal, signal_allcmpts, time, itertimes, totaltime)
end
