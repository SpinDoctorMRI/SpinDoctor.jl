# Bloch-Torrey PDE solvers
#
# The semi-discrete Bloch-Torrey equation is M u' = -J(t) u with J(t) = A + i gamma g(t) . Mx,
# where A = S + Q + R (stiffness, flux and relaxation matrices), M is the mass matrix and Mx
# the first order moment matrices. Further time stepping methods (OrdinaryDiffEq algorithms)
# are available in the SciMLBase extension.

"""
    ThetaSolver(; theta = 0.5, timestep = nothing)

Theta-scheme time stepping for the BTPDE (`theta = 0.5`: Crank-Nicolson, second order;
`theta = 1`: implicit Euler, first order), with steps of at most `timestep` in each
characteristic interval of the gradient (by default a thousandth of the echo time). The
gradient is evaluated at the middle of each step.

For interval-wise constant scalar gradients (e.g. [`PGSE`](@ref), [`DoublePGSE`](@ref)),
the linear systems are solved with one sparse factorization per interval. For other
gradients, they are solved with GMRES, preconditioned with a factorization of the diffusion
part (one per step size), so that the matrices are not refactorized at every step.
"""
Base.@kwdef struct ThetaSolver
    theta::Float64 = 0.5
    timestep::Union{Nothing,Float64} = nothing
end

"""
    solve_btpde(model, matrices, gradient, solver = ThetaSolver(); callbacks = [])

Solve the Bloch-Torrey partial differential equation using P1 finite elements in space and
the time stepping `solver`, and return the magnetization at the echo time.

The time stepping is a [`ThetaSolver`](@ref). With OrdinaryDiffEq, its ODE algorithms may
also be used, e.g. `using OrdinaryDiffEqBDF; solve_btpde(model, matrices, gradient, QNDF();
abstol, reltol)`.

The `callbacks` (see [`AbstractCallback`](@ref)) are called after each time step.
"""
function solve_btpde(
    model::Model{T,dim},
    matrices,
    gradient,
    solver::ThetaSolver = ThetaSolver();
    callbacks = AbstractCallback[],
) where {T,dim}
    problem = (; model, matrices)
    (; theta) = solver
    (; gamma) = model
    (; M, S, R, Mx, Q) = matrices
    timestep = something(solver.timestep, echotime(gradient) / 1000)
    A = S + Q + R
    J(t) = A + im * gamma * sum(gradient(t)[d] * Mx[d] for d = 1:dim)
    constant = gradient isa ScalarGradient && isconstant(gradient.profile)

    t = zero(T)
    u = initial_conditions(model)
    for cb in callbacks
        initialize!(cb, problem, gradient, u, t)
    end

    # Preconditioners: factorizations of the diffusion part, one per step size
    preconditioners = Dict{T,Any}()
    preconditioner(h) = get!(preconditioners, h) do
        P = M + h * theta * A
        issymmetric(P) ? cholesky(Symmetric(P)) : lu(P)
    end

    ivals = intervals(gradient)
    for i = 1:(length(ivals)-1)
        nstep = max(1, ceil(Int, (ivals[i+1] - ivals[i]) / timestep))
        h = (ivals[i+1] - ivals[i]) / nstep
        if constant
            # Constant gradient in the interval: one factorization
            Ji = J((ivals[i] + ivals[i+1]) / 2)
            F = lu(M + h * theta * Ji)
            E = M - h * (1 - theta) * Ji
        else
            P = preconditioner(h)
        end
        for n = 1:nstep
            t0 = ivals[i] + (n - 1) * h
            t1 = n == nstep ? ivals[i+1] : t0 + h
            if constant
                u = F \ (E * u)
            else
                Jn = J((t0 + t1) / 2)
                rhs = M * u - h * (1 - theta) * (Jn * u)
                u = preconditioned_gmres(M + h * theta * Jn, rhs, u, P)
            end
            t = t1
            for cb in callbacks
                update!(cb, problem, gradient, u, t)
            end
        end
    end

    finalize!.(callbacks)

    u
end

# Solve K x = b with GMRES, left preconditioned with the real factorization P
function preconditioned_gmres(K, b, x0, P)
    Pinv(v) = (w = P \ [real(v) imag(v)]; complex.(w[:, 1], w[:, 2]))
    x, info = KrylovKit.linsolve(
        v -> Pinv(K * v),
        Pinv(b),
        x0,
        KrylovKit.GMRES(; krylovdim = 30, maxiter = 100, tol = 1e-10 * norm(x0)),
    )
    info.converged > 0 || @warn "GMRES did not converge" info.normres
    x
end
