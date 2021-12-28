"""
    run(simulation::IntervalConstanBTPDE, gradient)

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
This function uses a manual time stepping scheme (theta-rule), that requires a degree of
implicitness `θ` and a time step `Δt`.
    `θ = 0.5`: Crank-Nicolson (second order)
    `θ = 1.0`: Implicit Euler (first order)
"""
function solve(
    simulation::IntervalConstanBTPDE{T},
    gradient::ScalarGradient;
    callbacks = AbstractCallback[],
) where {T}
    @unpack θ, timestep, model, matrices = simulation
    @unpack mesh, D, T₂, ρ, γ = model
    @unpack M, S, R, Mx, Q, M_cmpts = matrices

    isconstant(gradient.profile) || error("Time profile must be interval-wise constant")
    ivals = intervals(gradient.profile)

    # Create initial conditions (enforce complex values)
    ρ = mapreduce((ρ, p) -> fill(ρ, size(p, 2)), vcat, ρ, mesh.points)

    Jac!(J, g⃗) =
           @.(J = -(S + Q + R + im * γ * (g⃗[1] * Mx[1] + g⃗[2] * Mx[2] + g⃗[3] * Mx[3])))

    t = 0.0
    ξ = copy(ρ)
    Ey = copy(ξ)
    J = -(S + Q + im * sum(Mx)) # Allocate sparsity pattern

    for cb ∈ callbacks
        initialize!(cb, simulation, gradient, ξ, t)
    end

    # Crank-Nicolson time stepping
    for i = 1:(length(ivals) - 1)
        @debug "Solving for interval [%g, %g]\n" ivals[i] ivals[i + 1]

        # Adjust time step to divide interval uniformly
        ival_length = ivals[i + 1] - ivals[i]
        nt = round(Int, ival_length / timestep)
        Δt = ival_length / nt

        # Compute gradient at midpoint of interval
        tmid = (ivals[i] + ivals[i + 1]) / 2
        g⃗ = gradient(tmid)

        # Build matrices for interval
        Jac!(J, g⃗)
        F = lu(M .- Δt .* θ .* J)
        # F = factorize(M .- Δt .* θ .* J)
        E = @. complex(M) + Δt * (1 - θ) * J

        # Step to end of interval
        for it = 1:nt
            mul!(Ey, E, ξ)
            ldiv!(ξ, F, Ey)
            t += Δt
            for cb ∈ callbacks
                update!(cb, simulation, gradient, ξ, t)
            end
        end
    end

    finalize!.(callbacks)

    ξ
end