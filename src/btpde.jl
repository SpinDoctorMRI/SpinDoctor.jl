"""
    IntervalConstantSolver(; θ = 0.5, timestep)

BTPDE solver specialized on intervalwise constant `ScalarGradient`s, e.g [`PGSE`](@ref),
[`DoublePGSE`](@ref). Raises an error for other gradients.
"""
Base.@kwdef struct IntervalConstantSolver{T}
    θ::T = 0.5
    timestep::T
end

"""
    solve_btpde(
        model,
        matrices,
        gradient,
        odesolver = QNDF();
        abstol = 1e-6,
        reltol = 1e-4,
        callbacks = [],
    )

Solve the Bloch-Torrey partial differential equation using P1 finite elements in space and
the ODE algorithm `odesolver` in time. Return the magnetization at the echo time.

For interval-wise constant gradients, an [`IntervalConstantSolver`](@ref) may be used
instead of an ODE algorithm.
"""
function solve_btpde(
    model::Model{T,dim},
    matrices,
    gradient,
    odesolver::AbstractODEAlgorithm = QNDF();
    abstol = 1e-6,
    reltol = 1e-4,
    callbacks = AbstractCallback[],
) where {T,dim}
    problem = (; model, matrices)
    (; γ) = model
    (; M, S, R, Mx, Q) = matrices

    ρ = initial_conditions(model)

    # Time dependent Jacobian of ODE function with respect to the state `ξ`
    function jac!(J, ξ, p, t)
        (; S, Q, R, gradient) = p
        # println(t)
        g = gradient(t)
        if dim == 2
            @. J = -(S + Q + R + im * γ * (g[1] * Mx[1] + g[2] * Mx[2]))
        elseif dim == 3
            @. J = -(S + Q + R + im * γ * (g[1] * Mx[1] + g[2] * Mx[2] + g[3] * Mx[3]))
        end
    end

    # Time dependent ODE function
    function Mdξ!(dξ, ξ, p, t)
        jac!(p.J, ξ, p, t)
        mul!(dξ, p.J, ξ)
    end

    function func(u, t, integrator)
        for cb ∈ callbacks
            update!(cb, problem, gradient, u, t)
        end
    end

    # Jacobian sparsity pattern
    jac_prototype = -(S + Q + im * sum(Mx))

    # ODE problem
    J = copy(jac_prototype)
    p = (; J, S, Q, R, gradient)
    TE = echotime(gradient)

    for cb ∈ callbacks
        initialize!(cb, problem, gradient, ρ, 0.0)
    end

    odefunction = ODEFunction(Mdξ!; jac = jac!, mass_matrix = M, jac_prototype)
    odeproblem = ODEProblem(odefunction, ρ, (0, TE), p)

    callback = FunctionCallingCallback(func; func_start = false)

    # Solve ODE problem
    sol = SciMLBase.solve(
        odeproblem,
        odesolver;
        reltol,
        abstol,
        callback,
        save_everystep = false,
        tstops = intervals(gradient)[2:end],
    )

    finalize!.(callbacks)

    sol.u[end]
end

"""
    solve_btpde(
        model,
        matrices,
        gradient::ScalarGradient,
        solver::IntervalConstantSolver;
        callbacks = [],
    )

Solve the Bloch-Torrey partial differential equation using P1 finite elements in space and a
theta-rule in time.
This time stepping scheme requires a degree of
implicitness `θ` and a time step `Δt`:

  - `θ = 0.5`: Crank-Nicolson (second order)
  - `θ = 1.0`: Implicit Euler (first order)

The function only works for interval-wise constant `ScalarGradient`s, and errors otherwise.
"""
function solve_btpde(
    model::Model{T,dim},
    matrices,
    gradient::ScalarGradient,
    solver::IntervalConstantSolver{T};
    callbacks = AbstractCallback[],
) where {T,dim}
    problem = (; model, matrices)
    (; θ, timestep) = solver
    (; γ) = model
    (; M, S, R, Mx, Q) = matrices

    isconstant(gradient.profile) || error("Time profile must be interval-wise constant")
    ivals = intervals(gradient.profile)

    ρ = initial_conditions(model)

    function jac!(J, g)
        if dim == 2
            @.(J = -(S + Q + R + im * γ * (g[1] * Mx[1] + g[2] * Mx[2])))
        elseif dim == 3
            @.(J = -(S + Q + R + im * γ * (g[1] * Mx[1] + g[2] * Mx[2] + g[3] * Mx[3])))
        end
    end

    t = 0.0
    ξ = copy(ρ)
    Ey = copy(ξ)
    J = -(S + Q + im * sum(Mx)) # Allocate sparsity pattern

    for cb ∈ callbacks
        initialize!(cb, problem, gradient, ξ, t)
    end

    # Crank-Nicolson time stepping
    for i = 1:(length(ivals)-1)
        @debug "Solving for interval [%g, %g]\n" ivals[i] ivals[i+1]

        # Adjust time step to divide interval uniformly
        ival_length = ivals[i+1] - ivals[i]
        nt = round(Int, ival_length / timestep)
        Δt = ival_length / nt

        # Compute gradient at midpoint of interval
        tmid = (ivals[i] + ivals[i+1]) / 2
        g = gradient(tmid)

        # Build matrices for interval
        jac!(J, g)
        F = lu(M .- Δt .* θ .* J)
        # F = factorize(M .- Δt .* θ .* J)
        E = @. complex(M) + Δt * (1 - θ) * J

        # Step to end of interval
        for _ = 1:nt
            mul!(Ey, E, ξ)
            ldiv!(ξ, F, Ey)
            t += Δt
            for cb ∈ callbacks
                update!(cb, problem, gradient, ξ, t)
            end
        end
    end

    finalize!.(callbacks)

    ξ
end
