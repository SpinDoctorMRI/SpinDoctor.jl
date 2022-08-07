"""
    solve(
        problem::HADC,
        gradient::ScalarGradient,
        odesolver = QNDF();
        abstol = 1e-6,
        reltol = 1e-4,
    )

Compute the ADC using a homogenized ADC model (HADC). This is currently only
implemented for scalar gradients.
"""
function solve(
    problem::HADC,
    gradient::ScalarGradient,
    odesolver = QNDF(; autodiff = false);
    abstol = 1e-6,
    reltol = 1e-4,
)
    (; model, matrices) = problem
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
        # @show t
        (; mS, f, surfint) = p
        mul!(dω, mS, ω)
        dω .+= integral(f, t) .* surfint
    end

    # Create ODE function and Jacobian from matrices
    Jac!(J, _, p, t) = (J .= p.mS)

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
            jac = Jac!,
            jac_prototype = p.mS,
            mass_matrix = M_cmpts[icmpt],
        )
        odeproblem = ODEProblem(odefunction, ω₀[icmpt], (0, TE), p; progress = false)

        # Solve ODE, keep all time steps (for integral)
        sol = OrdinaryDiffEq.solve(
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
