"""
    solve(problem::GeneralBTPDE, gradient[; callbacks])

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
"""
function solve(
    problem::GeneralBTPDE{T},
    gradient;
    callbacks = AbstractCallback[],
) where {T}
    (; model, matrices, reltol, abstol, odesolver) = problem
    (; mesh, D, T₂, ρ, γ) = model
    (; M, S, R, Mx, Q, M_cmpts) = matrices

    ρ = initial_conditions(model)

    # Time dependent ODE function
    function Mdξ!(dξ, ξ, p, t)
        Jac!(p.J, ξ, p, t)
        mul!(dξ, p.J, ξ)
    end

    # Time dependent Jacobian of ODE function with respect to the state `ξ`
    function Jac!(J, ξ, p, t)
        (; S, Q, R, gradient) = p
        g⃗ = gradient(t)
        @. J = -(S + Q + R + im * γ * (g⃗[1] * Mx[1] + g⃗[2] * Mx[2] + g⃗[3] * Mx[3]))
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
        initialize!(cb, problem, gradient, ρ, 0)
    end

    odefunction = ODEFunction(Mdξ!; mass_matrix = M, jac = Jac!, jac_prototype)
    odeproblem = ODEProblem(odefunction, ρ, (0, TE), p, progress = false)

    callback = FunctionCallingCallback(func; func_start = false)

    # Solve ODE problem
    sol = OrdinaryDiffEq.solve(
        odeproblem,
        odesolver;
        reltol,
        abstol,
        callback,
        save_everystep = false,
    )

    finalize!.(callbacks)

    sol.u[end]
end
