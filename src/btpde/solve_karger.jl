"""
    solve_karger(model, difftensors, gradient; timestep)

Solve the finite pulse Karger model (FPK) using precomputed effective diffusion tensors
`difftensors` (one per compartment), and return the compartment signals.

The compartment signals are integrated with a fourth order Magnus method, with steps of at
most `timestep` within each characteristic interval of the gradient sequence.
"""
function solve_karger(model, difftensors, gradient::ScalarGradient; timestep)
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

    ADC_diag = Diagonal([d' * D * d for D ∈ difftensors])

    # Linear ODE dS/dt = J(t) S with J(t) = A - (γ g F(t))² ADC - R
    A = Matrix(A)
    R = Matrix(R)
    J(t) = A - (γ * g * integral(f, t))^2 * ADC_diag - R

    # Fourth order Magnus method (two Gauss-Legendre points per step)
    c₁, c₂ = 1 / 2 - √3 / 6, 1 / 2 + √3 / 6
    S = complex(S₀)
    ivals = intervals(gradient)
    for i = 1:(length(ivals)-1)
        nstep = max(1, ceil(Int, (ivals[i+1] - ivals[i]) / dt))
        h = (ivals[i+1] - ivals[i]) / nstep
        for istep = 1:nstep
            t = ivals[i] + (istep - 1) * h
            J₁ = J(t + c₁ * h)
            J₂ = J(t + c₂ * h)
            Ω = h / 2 * (J₁ + J₂) + √3 / 12 * h^2 * (J₂ * J₁ - J₁ * J₂)
            S = exp(Ω) * S
        end
    end

    S
end
