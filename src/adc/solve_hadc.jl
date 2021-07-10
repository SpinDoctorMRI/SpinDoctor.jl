"""
    solve_hadc(model, experiment) -> NamedTuple

Compute ADCs using the homogenized ADC model (HADC).
"""
function solve_hadc(model, experiment)

    # Measure function evalutation time
    starttime = Base.time()

    # Extract input parameters
    @unpack mesh, D, T₂, ρ = model
    @unpack directions, sequences = experiment.gradient
    @unpack odesolver, reltol, abstol = experiment.hadc

    # Deduce sizes
    ncompartment, nboundary = size(mesh.facets)
    npoint_cmpts = size.(mesh.points, 2)
    inds_cmpts = cumsum([0; npoint_cmpts])
    ndirection = size(directions, 2)
    nsequence = length(sequences)

    # Number of points in each compartment
    npoint_cmpts = size.(mesh.points, 2)

    # Assemble finite element matrices
    M = []
    S = []
    G = []
    volumes = zeros(ncompartment)
    for icmpt = 1:ncompartment
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        fevolumes, = get_mesh_volumes(points, elements)
        volumes[icmpt] = sum(fevolumes)

        # Assemble finite element matrices
        push!(M, assemble_mass_matrix(elements', fevolumes))
        push!(S, assemble_stiffness_matrix(elements', points', D[icmpt]))

        # Compute surface integrals
        push!(G, zeros(npoint_cmpts[icmpt], 3))
        for iboundary = 1:nboundary
            facets = mesh.facets[icmpt, iboundary]
            if !isempty(facets)
                # Get facet normals
                _, _, normals = get_mesh_surfacenormals(points, elements, facets)

                # Surface normal weighted flux matrix (in each canonical direction)
                for dim = 1:3
                    Q = assemble_flux_matrix(facets', points', normals[dim, :])
                    G[icmpt][:, dim] += sum(Q, dims = 2)
                end
            end
        end
    end

    # Initial conditions
    ω₀ = zeros.(npoint_cmpts)

    # Time dependent ODE function
    function Mdω!(dω, ω, p, t)
        # @show t
        @unpack mS, f, surfint = p
        mul!(dω, mS, ω)
        dω .+= integral(f, t) .* surfint
    end

    # Allocate output arrays
    adc = zeros(ncompartment, nsequence, ndirection)
    adc_allcmpts = zeros(nsequence, ndirection)
    itertimes = zeros(ncompartment, nsequence, ndirection)

    # Iterate over compartments and gradient sequences and directions
    for idir = 1:ndirection, iseq = 1:nsequence, icmpt = 1:ncompartment

        # Measure iteration time
        itertime = time()

        # Extract parameters for iteration
        g = directions[:, idir]
        f = sequences[iseq]
        TE = echotime(f)

        # Free diffusivity in direction g
        D₀ = g' * D[icmpt] * g

        # Compute surface integrals in gradient direction
        surfint = G[icmpt] * (D[icmpt] * g)

        # Display state of iterations
        println("Solving HADC model of size $(sum(npoint_cmpts)):")
        println("  Direction   $idir of $ndirection: g = $g")
        println("  Sequence    $iseq of $nsequence: f = $f")
        println("  Compartment $icmpt of $ncompartment")

        # Create ODE function and Jacobian from matrices
        jac = (J, _, p, t) -> (J .= p.mS)
        p = (; mS = -S[icmpt], f, surfint)
        odefunction =
            ODEFunction(Mdω!, jac = jac, jac_prototype = p.mS, mass_matrix = M[icmpt])
        odeproblem = ODEProblem(odefunction, ω₀[icmpt], (0, TE), p, progress = false)

        # Solve ODE, keep all time steps (for integral)
        sol = solve(odeproblem, odesolver, reltol = reltol, abstol = abstol)

        # Integral over compartment boundary
        a, err = quadgk(t -> integral(f, t) * (surfint' * sol(t)), 0, TE)

        # HADC (free diffusivity minus correction)
        adc[icmpt, iseq, idir] = D₀ - a / volumes[icmpt] / bvalue_no_q(f)

        # Computational time
        itertimes[icmpt, iseq, idir] = time() - itertime
    end

    # Compute total HADC (weighted sum over compartments)
    weights = ρ .* volumes
    weights = weights / sum(weights)
    adc_allcmpts[:] = sum(weights .* adc, dims = 1)

    # Create output structure
    totaltime = time() - starttime

    (; adc, adc_allcmpts, totaltime, itertimes)
end
