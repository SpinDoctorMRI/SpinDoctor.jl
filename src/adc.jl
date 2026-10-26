"""
    compute_adc_sta(model)

Compute the ADC for each compartment in the short diffusion time regime.
"""
function compute_adc_sta(model::Model{T}, gradient::ScalarGradient) where {T}
    ncompartment = length(model.mesh.points)

    adcs = zeros(T, ncompartment)
    for icmpt = 1:ncompartment
        points = model.mesh.points[icmpt]
        elements = model.mesh.elements[icmpt]
        facets = reduce(hcat, model.mesh.facets[icmpt, :])

        volumes, = get_mesh_volumes(points, elements)
        v = sum(volumes)
        areas, _, normals = get_mesh_surface(points, facets)
        d = gradient.dir
        D₀ = d' * model.D[icmpt] * d

        # Project surface areas onto plane orthogonal to gradient direction
        # (result is 1×1 matrix, extract scalar component)
        SAu = ((d'normals) .^ 2*areas)[1]

        t = diffusion_time_sta(gradient.profile)
        adcs[icmpt] = D₀ * (1 - 4 / 3 * √(D₀ / π) * SAu / v * √t)
    end

    all(≥(0), adcs) ||
        @warn "Obtained negative STA ADC, STA approximation does not hold" filter(
            <(0),
            adcs,
        )

    adcs
end

"""
    solve_hadc(model, matrices, gradient, odesolver = QNDF(); abstol = 1e-6, reltol = 1e-4)

Compute the ADC of each compartment using the homogenized ADC model (HADC). This is
currently only implemented for scalar gradients.

The HADC model assumes impermeable compartment boundaries: permeabilities, initial spin
densities and T₂ relaxation are not taken into account. The ADC of a sample may be
estimated by weighting the compartment ADCs with [`compartment_weights`](@ref).
"""
function solve_hadc(
    model,
    matrices,
    gradient::ScalarGradient,
    odesolver = QNDF();
    abstol = 1e-6,
    reltol = 1e-4,
)
    (; mesh, D) = model
    (; M_cmpts, S_cmpts, G, volumes) = matrices

    f = gradient.profile
    dir = gradient.dir
    TE = echotime(f)

    # Deduce sizes
    ncompartment = length(D)

    # Initial conditions
    ω₀ = zeros.(size.(mesh.points, 2))

    # Time dependent ODE function
    function Mdω!(dω, ω, p, t)
        (; mS, f, surfint) = p
        mul!(dω, mS, ω)
        dω .+= integral(f, t) .* surfint
    end

    # Create ODE function and Jacobian from matrices
    jac!(J, ω, p, t) = (J .= p.mS)

    # Allocate output array
    adc = zeros(ncompartment)

    # Iterate over compartments and gradient sequences and directions
    for icmpt = 1:ncompartment
        # Free diffusivity in gradient direction
        D₀ = dir' * D[icmpt] * dir

        # Surface integrals in gradient direction
        surfint = G[icmpt] * (D[icmpt] * dir)

        p = (; mS = -S_cmpts[icmpt], f, surfint)

        odefunction = ODEFunction(
            Mdω!;
            jac = jac!,
            jac_prototype = p.mS,
            mass_matrix = M_cmpts[icmpt],
        )
        odeproblem = ODEProblem(odefunction, ω₀[icmpt], (0, TE), p)

        # Solve ODE, keep all time steps (for integral)
        sol = SciMLBase.solve(
            odeproblem,
            odesolver;
            reltol = reltol,
            abstol = abstol,
            tstops = intervals(gradient)[2:end],
        )

        # Integral over compartment boundary
        a, = quadgk(t -> integral(f, t) * (surfint' * sol(t)), 0, TE)

        # HADC (free diffusivity minus correction)
        adc[icmpt] = D₀ - a / volumes[icmpt] / int_F²(f)
    end

    adc
end
