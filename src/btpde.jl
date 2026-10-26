# Bloch-Torrey PDE solvers

"""
    IntervalConstantSolver(; theta = 0.5, timestep)

Time stepping for the BTPDE with interval-wise constant scalar gradients, e.g.
[`PGSE`](@ref) and [`DoublePGSE`](@ref): a theta-scheme with steps of at most `timestep`
in each interval (`theta = 0.5`: Crank-Nicolson, `theta = 1`: implicit Euler). See
[`solve_btpde`](@ref).
"""
Base.@kwdef struct IntervalConstantSolver{T}
    theta::T = 0.5
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
    (; gamma) = model
    (; M, S, R, Mx, Q) = matrices

    xi0 = initial_conditions(model)

    # Time dependent Jacobian of ODE function with respect to the state `xi`
    function jac!(J, xi, p, t)
        (; S, Q, R, gradient) = p
        g = gradient(t)
        if dim == 2
            @. J = -(S + Q + R + im * gamma * (g[1] * Mx[1] + g[2] * Mx[2]))
        elseif dim == 3
            @. J = -(S + Q + R + im * gamma * (g[1] * Mx[1] + g[2] * Mx[2] + g[3] * Mx[3]))
        end
    end

    # Time dependent ODE function
    function Mdxi!(dxi, xi, p, t)
        jac!(p.J, xi, p, t)
        mul!(dxi, p.J, xi)
    end

    function func(u, t, integrator)
        for cb in callbacks
            update!(cb, problem, gradient, u, t)
        end
    end

    # Jacobian sparsity pattern
    jac_prototype = -(S + Q + im * sum(Mx))

    # ODE problem
    J = copy(jac_prototype)
    p = (; J, S, Q, R, gradient)
    TE = echotime(gradient)

    for cb in callbacks
        initialize!(cb, problem, gradient, xi0, zero(T))
    end

    odefunction = ODEFunction(Mdxi!; jac = jac!, mass_matrix = M, jac_prototype)
    odeproblem = ODEProblem(odefunction, xi0, (0, TE), p)

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

Solve the Bloch-Torrey partial differential equation using P1 finite elements in space and
the theta-scheme of `solver` in time. The time profile of the gradient must be interval-wise
constant.
"""
function solve_btpde(
    model::Model{T,dim},
    matrices,
    gradient::ScalarGradient,
    solver::IntervalConstantSolver{T};
    callbacks = AbstractCallback[],
) where {T,dim}
    problem = (; model, matrices)
    (; theta, timestep) = solver
    (; gamma) = model
    (; M, S, R, Mx, Q) = matrices

    isconstant(gradient.profile) || error("Time profile must be interval-wise constant")
    ivals = intervals(gradient.profile)

    function jac!(J, g)
        if dim == 2
            @.(J = -(S + Q + R + im * gamma * (g[1] * Mx[1] + g[2] * Mx[2])))
        elseif dim == 3
            @.(J = -(S + Q + R + im * gamma * (g[1] * Mx[1] + g[2] * Mx[2] + g[3] * Mx[3])))
        end
    end

    t = zero(T)
    xi = initial_conditions(model)
    Ey = copy(xi)
    J = -(S + Q + im * sum(Mx)) # Allocate sparsity pattern

    for cb in callbacks
        initialize!(cb, problem, gradient, xi, t)
    end

    # Theta-scheme time stepping
    for i = 1:(length(ivals)-1)
        # Divide the interval into steps of at most `timestep`
        ival_length = ivals[i+1] - ivals[i]
        nt = max(1, ceil(Int, ival_length / timestep))
        dt = ival_length / nt

        # Compute gradient at midpoint of interval
        tmid = (ivals[i] + ivals[i+1]) / 2
        g = gradient(tmid)

        # Build matrices for interval
        jac!(J, g)
        F = lu(M .- dt .* theta .* J)
        E = @. complex(M) + dt * (1 - theta) * J

        # Step to end of interval
        for _ = 1:nt
            mul!(Ey, E, xi)
            ldiv!(xi, F, Ey)
            t += dt
            for cb in callbacks
                update!(cb, problem, gradient, xi, t)
            end
        end
    end

    finalize!.(callbacks)

    xi
end
