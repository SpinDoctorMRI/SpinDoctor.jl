# Apparent diffusion coefficients (ADC)

"""
    compute_adc_sta(model, gradient)

Compute the ADC of each compartment in the short diffusion time approximation (STA), for
the direction of the scalar gradient `gradient`.
"""
function compute_adc_sta(model::Model{T}, gradient::ScalarGradient) where {T}
    (; mesh, D) = model
    d = gradient.dir
    t = diffusion_time_sta(gradient.profile)

    adcs = zeros(T, length(mesh.points))
    for icmpt in eachindex(mesh.points)
        points = mesh.points[icmpt]
        facets = reduce(hcat, mesh.facets[icmpt, :])
        volume = sum(first(element_geometry(points, mesh.elements[icmpt])))
        areas, _, normals = facet_geometry(points, facets)
        D0 = d' * D[icmpt] * d

        # Surface area projected onto the plane orthogonal to the gradient direction
        area = sum((d' * normals)' .^ 2 .* areas)

        adcs[icmpt] = D0 * (1 - 4 / 3 * sqrt(D0 / pi) * area / volume * sqrt(t))
    end

    all(>=(0), adcs) ||
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
densities and T2 relaxation are not taken into account. The ADC of a sample may be
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
    omega0 = zeros.(size.(mesh.points, 2))

    # Time dependent ODE function
    function Mdomega!(domega, omega, p, t)
        (; mS, f, surfint) = p
        mul!(domega, mS, omega)
        domega .+= integral(f, t) .* surfint
    end

    # Create ODE function and Jacobian from matrices
    jac!(J, omega, p, t) = (J .= p.mS)

    # Allocate output array
    adc = zeros(eltype(volumes), ncompartment)

    # Iterate over compartments and gradient sequences and directions
    for icmpt = 1:ncompartment
        # Free diffusivity in gradient direction
        D0 = dir' * D[icmpt] * dir

        # Surface integrals in gradient direction
        surfint = G[icmpt] * (D[icmpt] * dir)

        p = (; mS = -S_cmpts[icmpt], f, surfint)

        odefunction = ODEFunction(
            Mdomega!;
            jac = jac!,
            jac_prototype = p.mS,
            mass_matrix = M_cmpts[icmpt],
        )
        odeproblem = ODEProblem(odefunction, omega0[icmpt], (0, TE), p)

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
        adc[icmpt] = D0 - a / volumes[icmpt] / int_F2(f)
    end

    adc
end
