"""
    solve(problem::Karger, gradient)

Solve the finite pulse Karger model (FPK) using precomputed effective diffusion tensors
`difftensors`.
"""
function solve(problem::Karger, gradient::ScalarGradient)
    (; model, difftensors, odesolver, timestep) = problem
    (; mesh, T₂, ρ, κ, γ) = model

    dt = timestep
    f = gradient.profile
    d = gradient.dir
    g = gradient.amplitude
    TE = echotime(f)

    # Deduce sizes
    ncompartment, nboundary = size(mesh.facets)

    # Volumes
    volumes = get_cmpt_volumes(mesh)

    # Compute surface areas
    surface_areas = spzeros(ncompartment, nboundary)
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt]
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

    # Update function for linear ODE operator
    function update_func(J, u, p, t)
        (; A, ADC_diag, R, f, γ, g) = p
        J .= A .- (integral(f, t)^2 * (γ * g)^2) .* ADC_diag .- R
    end

    ADC_diag = spdiagm([d' * D * d for D ∈ difftensors])

    p = (; A, ADC_diag, R, f, γ, g)
    J_prototype = A - ADC_diag - R
    J = DiffEqArrayOperator(J_prototype; update_func)

    prob = ODEProblem(J, S₀, (0.0, TE), p)
    sol = OrdinaryDiffEq.solve(prob, odesolver; dt)

    sol.u[end]
end
