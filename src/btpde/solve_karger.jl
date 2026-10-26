"""
    solve(problem::Karger, gradient, odesolver = QNDF(); timestep)

Solve the finite pulse Karger model (FPK) using precomputed effective diffusion tensors
`difftensors`.
"""
function solve(problem::Karger, gradient::ScalarGradient, odesolver = QNDF(); timestep)
    (; model, difftensors) = problem
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
                areas = first(get_mesh_surface(points, facets[iboundary]))
                surface_areas[icmpt, iboundary] = sum(areas)
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

    ADC_diag = spdiagm([d' * D * d for D ∈ difftensors])

    # Time-dependent linear ODE: du/dt = J(t) * u
    # where J(t) = A - (F(t)² * (γg)²) * ADC_diag - R
    J_cache = copy(A - ADC_diag - R)
    p = (; A, ADC_diag, R, f, γ, g, J_cache)

    function karger_rhs!(du, u, p, t)
        (; A, ADC_diag, R, f, γ, g, J_cache) = p
        J_cache .= A .- (integral(f, t)^2 * (γ * g)^2) .* ADC_diag .- R
        mul!(du, J_cache, u)
    end

    prob = ODEProblem(karger_rhs!, S₀, (0.0, TE), p)
    sol = OrdinaryDiffEqCore.solve(prob, QNDF(); dt)

    sol.u[end]
end
