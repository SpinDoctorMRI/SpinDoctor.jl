"""
Solve the BTPDE with the ODE algorithms of OrdinaryDiffEq (loaded with an OrdinaryDiffEq
solver package, e.g. `using OrdinaryDiffEqBDF` for `QNDF`).
"""
module SpinDoctorSciMLExt

using LinearAlgebra
using SciMLBase: SciMLBase, AbstractODEAlgorithm, DiscreteCallback, ODEFunction, ODEProblem
using SpinDoctor
using SpinDoctor:
    AbstractCallback, Model, echotime, finalize!, initialize!, intervals, update!
import SpinDoctor: solve_btpde

function solve_btpde(
    model::Model{T,dim},
    matrices,
    gradient,
    odesolver::AbstractODEAlgorithm;
    abstol = 1e-6,
    reltol = 1e-4,
    callbacks = AbstractCallback[],
) where {T,dim}
    problem = (; model, matrices)
    (; gamma) = model
    (; M, S, R, Mx, Q) = matrices
    A = S + Q + R

    # Time dependent Jacobian of the ODE function
    function jac!(J, u, p, t)
        g = gradient(t)
        if dim == 2
            @. J = -(A + im * gamma * (g[1] * Mx[1] + g[2] * Mx[2]))
        else
            @. J = -(A + im * gamma * (g[1] * Mx[1] + g[2] * Mx[2] + g[3] * Mx[3]))
        end
    end

    # ODE function M u' = J(t) u
    function f!(du, u, p, t)
        jac!(p.J, u, p, t)
        mul!(du, p.J, u)
    end

    # Jacobian sparsity pattern
    jac_prototype = -(A + im * sum(Mx))
    p = (; J = copy(jac_prototype))

    u0 = initial_conditions(model)
    for cb in callbacks
        initialize!(cb, problem, gradient, u0, zero(T))
    end
    function affect!(integrator)
        for cb in callbacks
            update!(cb, problem, gradient, integrator.u, integrator.t)
        end
        SciMLBase.u_modified!(integrator, false)
    end
    callback = DiscreteCallback(
        (u, t, integrator) -> true,
        affect!;
        save_positions = (false, false),
    )

    odefunction = ODEFunction(f!; jac = jac!, mass_matrix = M, jac_prototype)
    odeproblem = ODEProblem(odefunction, u0, (zero(T), echotime(gradient)), p)
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

end
