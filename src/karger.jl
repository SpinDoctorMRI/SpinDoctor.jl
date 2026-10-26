# Finite pulse Karger model

"""
    solve_karger(model, difftensors, gradient; timestep)

Solve the finite pulse Karger model (FPK) using precomputed effective diffusion tensors
`difftensors` (one per compartment), and return the compartment signals.

The compartments exchange magnetization through their interfaces with rates given by the
permeabilities, interface areas and volumes. The compartment signals are integrated with a
fourth order Magnus method, with steps of at most `timestep` within each characteristic
interval of the gradient sequence.
"""
function solve_karger(model, difftensors, gradient::ScalarGradient; timestep)
    (; mesh, T2, rho, kappa, gamma) = model
    (; dir, profile, amplitude) = gradient
    ncompartment, nboundary = size(mesh.facets)
    volumes = compartment_volumes(mesh)

    # Exchange rates: an interface of area A and permeability kappa between compartments
    # i and j exchanges kappa A (m_j - m_i), with m = S / V the mean magnetization
    exchange = zeros(ncompartment, ncompartment)
    for iboundary = 1:nboundary
        cmpts = findall(!isempty, mesh.facets[:, iboundary])
        length(cmpts) == 2 || continue
        i, j = cmpts
        area = sum(first(facet_geometry(mesh.points[i], mesh.facets[i, iboundary])))
        exchange[i, j] = exchange[j, i] = kappa[iboundary] * area
    end
    A = exchange ./ volumes' - Diagonal(vec(sum(exchange; dims = 2)) ./ volumes)

    # Linear ODE dS/dt = J(t) S with J(t) = A - (gamma g F(t))^2 ADC - R
    adcs = Diagonal([dir' * D * dir for D in difftensors])
    R = Diagonal(1 ./ T2)
    J(t) = A - (gamma * amplitude * integral(profile, t))^2 * adcs - R

    # Fourth order Magnus method (two Gauss-Legendre points per step)
    c1, c2 = 1 / 2 - sqrt(3) / 6, 1 / 2 + sqrt(3) / 6
    S = complex(volumes .* rho)
    ivals = intervals(gradient)
    for i = 1:(length(ivals)-1)
        nstep = max(1, ceil(Int, (ivals[i+1] - ivals[i]) / timestep))
        h = (ivals[i+1] - ivals[i]) / nstep
        for istep = 1:nstep
            t = ivals[i] + (istep - 1) * h
            J1 = J(t + c1 * h)
            J2 = J(t + c2 * h)
            Omega = h / 2 * (J1 + J2) + sqrt(3) / 12 * h^2 * (J2 * J1 - J1 * J2)
            S = exp(Omega) * S
        end
    end

    S
end
