"""
    solve_karger(model, experiment, difftensors)
Solve the finite pulse Karger model (FPK) using precomputed effective diffusion tensors
`difftensors`.
"""
function solve_karger(model, experiment, difftensors)

    # Measure function evalutation time
    starttime = Base.time()

    # Extract input parameters
    @unpack mesh, T₂, ρ, κ = model
    @unpack directions, sequences, values, values_type = experiment.gradient
    @unpack odesolver, timestep = experiment.karger
    dt = timestep

    # Deduce sizes
    ncompartment, nboundary = size(mesh.facets)
    npoint_cmpts = size.(mesh.points, 2)
    inds_cmpts = cumsum([0; npoint_cmpts])
    ndirection = size(directions, 2)
    nsequence = length(sequences)
    namplitude = length(values)

    # Volumes
    volumes = get_cmpt_volumes(mesh)

    # Compute surface areas
    surface_areas = spzeros(ncompartment, nboundary)
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        facets = mesh.facets[icmpt, :]

        # Surface area
        for iboundary = 1:nboundary
            if !isempty(facets[iboundary])
                surface_areas[icmpt, iboundary], =
                    get_mesh_surface(points, facets[iboundary])
            end
        end
    end

    # Associate surface areas to the compartment interfaces they measure
    tau_inv = spzeros(ncompartment, ncompartment)
    for iboundary = 1:nboundary
        inds = surface_areas[:, iboundary].nzind
        if length(inds) == 2
            i₁ = inds[1]
            i₂ = inds[2]
            tmp = κ[iboundary] * surface_areas[i₁, iboundary]
            tau_inv[i₁, i₂] = tmp
            tau_inv[i₂, i₁] = tmp
        end
    end
    tau_inv = tau_inv ./ volumes'
    A = tau_inv - spdiagm((tau_inv * volumes) ./ volumes)

    # Relaxation tensor
    R = spdiagm(1 ./ T₂)

    # Initial signal
    S₀ = volumes .* ρ

    # Q-values and b-values
    if values_type == "q"
        qvalues = repeat(values, 1, nsequence)
        bvalues = values .^ 2 .* bvalue_no_q.(sequences)'
    else
        bvalues = repeat(values, 1, nsequence)
        qvalues = .√(values ./ bvalue_no_q.(sequences)')
    end

    # Update function for linear ODE operator
    function update_func(J, u, p, t)
        @unpack A, ADC_diag, R, f, q = p
        J .= A .- (integral(f, t)^2 * q^2) .* ADC_diag .- R
    end

    # Allocate output arrays
    signal = zeros(ncompartment, namplitude, nsequence, ndirection)
    signal_allcmpts = zeros(namplitude, nsequence, ndirection)
    itertimes = zeros(namplitude, nsequence, ndirection)

    # Iterate over gradient amplitudes, sequences and directions
    for idir = 1:ndirection, iseq = 1:nsequence, iamp = 1:namplitude

        # Measure iteration time
        itertime = Base.time()

        # Gradient amplitude
        q = qvalues[iamp, iseq]
        b = bvalues[iamp, iseq]

        # Time profile
        f = sequences[iseq]
        TE = echotime(f)

        # Gradient direction
        g = directions[:, idir]

        # Display state of iterations
        println("Computing Karger signal:")
        println("  Direction $idir of $ndirection: g = $g")
        println("  Sequence  $iseq of $nsequence: f = $f")
        println("  Amplitude $iamp of $namplitude: q = $q, b = $b")

        ADC_diag = spdiagm([g' * D * g for D ∈ difftensors[:, iseq]])

        p = (; A, ADC_diag, R, f, q)
        J_prototype = A - ADC_diag - R
        J = DiffEqArrayOperator(J_prototype, update_func = update_func)

        prob = ODEProblem(J, S₀, (0.0, TE), p)
        sol = solve(prob, odesolver, dt = dt)

        # Save final signal
        signal[:, iamp, iseq, idir] = sol.u[end]

        # Store timing
        itertimes[iamp, iseq, idir] = Base.time() - itertime
    end

    signal_allcmpts = sum(signal, dims = 1)[1, :, :, :]

    totaltime = Base.time() - starttime

    (; signal, signal_allcmpts, itertimes, totaltime)
end
