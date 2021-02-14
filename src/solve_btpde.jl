"""
    results = solve_btpde(mesh, domain, experiment, directions)

Solve the Bloch-Torrey partial differential equation.
"""
function solve_btpde(mesh, domain, experiment, directions)

    # Extract input parameters
    @unpack boundary_markers, compartments, boundaries, ncmpt, nboundary, diffusivity, relaxation, permeability, initial_density = domain
    @unpack ndir, flat_dirs, direction, sequences, values, values_type, btpde = experiment

    # Extract solver scheme
    @unpack odesolver, reltol, abstol, nsave = btpde

    # Deduce sizes
    nsequence = length(sequences)
    namplitude = length(values)
    npoint_cmpts = size.(mesh.points, 2)
    inds_cmpts = cumsum([0; npoint_cmpts])

    # Assemble finite element matrices compartment-wise
    fem_mat_cmpts = (
        M = [],
        S = [],
        Q = [],
        Mx = [[] for idim=1:3]
    )
    for icmpt = 1:ncmpt
        # Finite elements
        points = mesh.points[icmpt];
        facets = mesh.facets[icmpt, :];
        elements = mesh.elements[icmpt];
        volumes = get_mesh_volumes(points, elements);

        # Assemble mass, stiffness and flux matrices
        push!(fem_mat_cmpts.M, assemble_mass_matrix(elements', volumes))
        push!(fem_mat_cmpts.S, assemble_stiffness_matrix(elements', points', diffusivity[icmpt]))
        push!(fem_mat_cmpts.Q, assemble_flux_matrix_cmpt(points, facets, permeability))

        # Assemble first order product moment matrices
        for dim = 1:3
            push!(fem_mat_cmpts.Mx[dim], assemble_mass_matrix(elements', volumes, points[dim, :]))
        end
    end

    # Assemble global finite element matrices
    M = blockdiag(fem_mat_cmpts.M...)
    S = blockdiag(fem_mat_cmpts.S...)
    Q = couple_flux_matrix(mesh, fem_mat_cmpts.Q)
    Mx = [blockdiag(fem_mat_cmpts.Mx[dim]...) for dim = 1:3]

    # Create initial conditions (enforce complex values)
    ρ_cmpts = [fill(Complex(initial_density[icmpt]), npoint_cmpts[icmpt]) for icmpt = 1:ncmpt];

    # Initial spin density on entire domain
    ρ = vcat(ρ_cmpts...)

    # Allocate output arrays
    signal = fill(Array{ComplexF64}(undef, 0), ncmpt, namplitude, nsequence, ndir)
    signal_allcmpts = fill(Array{ComplexF64}(undef, 0), namplitude, nsequence, ndir)
    magnetization = fill(Array{ComplexF64, 2}(undef, 0, 0), ncmpt, namplitude, nsequence, ndir)
    time = fill(Float64[], namplitude, nsequence, ndir)

    # Q-values and b-values
    if values_type == 'q'
        qvalues = repeat(values, 1, nsequence)
        bvalues = values.^2 .* bvalue_no_q.(sequences')
    else
        bvalues = repeat(values, 1, nsequence)
        qvalues = .√(values ./ bvalue_no_q.(sequences'))
    end

    save_everystep = nsave > 1

    # ODE function
    function M∂u∂t!(du, u, p, t)
        J, S, Q, iA, q, f = p
        @. J = -S - Q - im * f(t) * q * A
        mul!(du, J, u)
        nothing
    end

    # Jacobian of ODE function with respect to the state `u`
    function Jac!(J, u, p, t)
        TMP, S, Q, A, q, f = p
        @. J = -(S + Q + im * f(t) * q * A)
        nothing
    end

    # Gather ODE function
    odefunction = ODEFunction(
        M∂u∂t!,
        jac = Jac!,
        jac_prototype = -complex(S + Q),
        mass_matrix = M,
    )


    # Iterate over gradient amplitudes, time profiles and directions
    for (iamp, iseq, idir) ∈ Iterators.product(1:namplitude, 1:nsequence, 1:ndir)

        # Gradient amplitude
        q = qvalues[iamp, iseq]
        b = bvalues[iamp, iseq]

        # Time profile
        f = sequences[iseq]
        interval = (0, echotime(f))
        if nsave == 1
            saveat = interval[2]
            saveat = [interval[1], interval[2]]
        else
            saveat = LinRange(interval..., nsave)
        end

        # Gradient direction
        dir = directions[:, idir]

        # Display state of iterations
        @printf "Solving BTPDE with size %d\n" (sum(npoint_cmpts))
        @printf "  Direction %d of %d: g = [%.2f, %.2f, %.2f]\n" idir ndir dir...
        @printf "  Sequence %d of %d: f = %s\n" iseq nsequence f
        @printf "  Amplitude %d of %d: q = %g, b = %g\n" iamp namplitude q b

        # Gradient direction dependent finite element matrix
        A = sum(dir[dim] * Mx[dim] for dim = 1:3)

        # ODE problem
        J = complex(S + Q + A)
        p = (J, S, Q, A, q, f)
        odeproblem = ODEProblem(
            odefunction, ρ, interval, p,
            progress = true,
            reltol = reltol,
            abstol = abstol,
        )

        # Solve ODE problem
        if isnothing(odesolver)
            sol = solve(odeproblem, saveat=saveat, tstops = [f.δ, f.Δ])
        else
            sol = solve(odeproblem, odesolver, saveat=saveat, tstops = [f.δ, f.Δ])
        end
        @show size(sol.u)

        # Extract solution
        time[iamp, iseq, idir] = sol.t
        mag = hcat(sol.u...)

        # Split solution into compartments
        for icmpt = 1:ncmpt
            inds = inds_cmpts[icmpt]+1:inds_cmpts[icmpt+1]

            # Store magnetization in compartment
            magnetization[icmpt, iamp, iseq, idir] = mag[inds, :]

            # Integrate magnetization over compartment
            signal[icmpt, iamp, iseq, idir] = sum(fem_mat_cmpts.M[icmpt] * mag[inds, :]; dims=1)[:]

            # Compute total signal
            if icmpt == 1
                signal_allcmpts[iamp, iseq, idir] = signal[icmpt, iamp, iseq, idir]
            else
                signal_allcmpts[iamp, iseq, idir] .+= signal[icmpt, iamp, iseq, idir]
            end

        end # Compartments
    end # Iterations

    # Return named tuple
    (; magnetization, signal, signal_allcmpts, time)

end # Function
