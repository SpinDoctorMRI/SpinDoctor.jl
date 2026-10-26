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
