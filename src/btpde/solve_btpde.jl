"""
    solve(simulation::GeneralBTPDE, gradient[; callbacks])

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
"""
function solve(
    simulation::GeneralBTPDE{T},
    gradient;
    callbacks = AbstractCallback[],
) where {T}
    @unpack model, matrices, reltol, abstol, odesolver = simulation
    @unpack mesh, D, T₂, ρ, γ = model
    @unpack M, S, R, Mx, Q, M_cmpts = matrices

    # Create initial conditions (enforce complex values)
    ρ = mapreduce((ρ, p) -> fill(ρ, size(p, 2)), vcat, ρ, mesh.points)

    # Time dependent ODE function
    function Mdξ!(dξ, ξ, p, t)
        @show t
        Jac!(p.J, ξ, p, t)
        mul!(dξ, p.J, ξ)
    end

    # Time dependent Jacobian of ODE function with respect to the state `ξ`
    function Jac!(J, ξ, p, t)
        @unpack S, Q, R, gradient = p
         g⃗ = gradient(t)
           @. J = -(S + Q + R + im * γ * (g⃗[1] * Mx[1] + g⃗[2] * Mx[2] + g⃗[3] * Mx[3]))
    end

    function func(u, t, integrator)
        for cb ∈ callbacks
            update!(cb, simulation, gradient, u, t)
        end
    end

    # Jacobian sparsity pattern
    jac_prototype = -(S + Q + im * sum(Mx))

    # ODE problem
    J = copy(jac_prototype)
    p = (; J, S, Q, R, gradient)
    TE = echotime(gradient)

    for cb ∈ callbacks
        initialize!(cb, simulation, gradient, ρ, 0)
    end

    odefunction = ODEFunction(Mdξ!; mass_matrix = M, jac = Jac!, jac_prototype)
    odeproblem = ODEProblem(odefunction, ρ, (0, TE), p, progress = false)

    # tstops = intervals(gradient.profile)[2:(end - 1)]

    # t_affect = LinRange(0, TE, callbacks[1].naffect)[2:end]
    # callback = PresetTimeCallback(t_affect, affect!; save_positions = (false, false))

    # condition(u,t,integrator) = true
    # callback = DiscreteCallback(condition, affect!; save_positions = (false, false))

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