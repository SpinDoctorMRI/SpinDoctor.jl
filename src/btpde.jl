# Bloch-Torrey PDE solvers
#
# The semi-discrete Bloch-Torrey equation is M u' = -J(t) u with J(t) = A + i gamma g(t) . Mx,
# where A = S + Q + R (stiffness, flux and relaxation matrices), M is the mass matrix and Mx
# the first order moment matrices. Further time stepping methods (OrdinaryDiffEq algorithms)
# are available in the SciMLBase extension.

"""
    ThetaSolver(; theta = 0.5, timestep = nothing, arraytype = Array)

Theta-scheme time stepping for the BTPDE (`theta = 0.5`: Crank-Nicolson, second order;
`theta = 1`: implicit Euler, first order), with steps of at most `timestep` in each
characteristic interval of the gradient (by default a thousandth of the echo time). The
gradient is evaluated at the middle of each step.

For interval-wise constant scalar gradients (e.g. [`PGSE`](@ref), [`DoublePGSE`](@ref)),
the linear systems are solved with one sparse factorization per interval. For other
gradients, they are solved with GMRES, preconditioned with a factorization of the diffusion
part (one per step size), so that the matrices are not refactorized at every step.

With `arraytype = CuArray` (after `using CUDA, CUDSS`), the time stepping runs on an NVIDIA
GPU, with the sparse direct solvers of cuDSS.
"""
Base.@kwdef struct ThetaSolver{A}
    theta::Float64 = 0.5
    timestep::Union{Nothing,Float64} = nothing
    arraytype::A = Array
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
    (; theta, arraytype) = solver
    (; gamma) = model
    (; M, S, R, Mx, Q) = matrices
    timestep = something(solver.timestep, echotime(gradient) / 1000)
    A = S + Q + R
    J(t) = A + im * gamma * sum(gradient(t)[d] * Mx[d] for d = 1:dim)
    constant = gradient isa ScalarGradient && isconstant(gradient.profile)

    # The callbacks get the magnetization on the CPU
    t = zero(T)
    u0 = initial_conditions(model)
    for cb in callbacks
        initialize!(cb, problem, gradient, u0, t)
    end
    u = to_device(arraytype, u0)

    # For gradients that are not constant in the intervals: the matrices on the device, and
    # preconditioners (factorizations of the diffusion part, one per step size)
    if !constant
        Md, Ad = to_device(arraytype, M), to_device(arraytype, A)
        Mxd = [to_device(arraytype, Mx[d]) for d = 1:dim]
    end
    preconditioners = Dict{T,Any}()
    preconditioner(h) = get!(preconditioners, h) do
        P = M + h * theta * A
        issymmetric(P) ? device_cholesky(arraytype, P) : device_lu(arraytype, P)
    end

    ivals = intervals(gradient)
    for i = 1:(length(ivals)-1)
        nstep = max(1, ceil(Int, (ivals[i+1] - ivals[i]) / timestep))
        h = (ivals[i+1] - ivals[i]) / nstep
        if constant
            # Constant gradient in the interval: one factorization
            Ji = J((ivals[i] + ivals[i+1]) / 2)
            F = device_lu(arraytype, M + h * theta * Ji)
            E = to_device(arraytype, M - h * (1 - theta) * Ji)
        else
            P = preconditioner(h)
        end
        for n = 1:nstep
            t0 = ivals[i] + (n - 1) * h
            t1 = n == nstep ? ivals[i+1] : t0 + h
            if constant
                u = F \ (E * u)
            else
                # J(t) changes at every step, so it is applied without being formed
                g = gradient((t0 + t1) / 2)
                Jmul(v) = Ad * v + im * gamma * sum(g[d] * (Mxd[d] * v) for d = 1:dim)
                rhs = Md * u - h * (1 - theta) * Jmul(u)
                u = preconditioned_gmres(v -> Md * v + h * theta * Jmul(v), rhs, u, P)
            end
            t = t1
            if !isempty(callbacks)
                ucpu = to_cpu(u)
                for cb in callbacks
                    update!(cb, problem, gradient, ucpu, t)
                end
            end
        end
    end

    finalize!.(callbacks)

    to_cpu(u)
end

# Arrays and sparse factorizations on the device of `arraytype` (the GPU methods are in the
# CUDA extension)
to_device(::Type{<:Array}, x) = x
device_lu(::Type{<:Array}, A) = lu(A)
device_cholesky(::Type{<:Array}, A) = cholesky(Symmetric(A))
to_device(arraytype, x) =
    error("Unsupported array type $arraytype (`CuArray` requires `using CUDA, CUDSS`)")
to_cpu(x) = x isa Array ? x : Array(x)

# Solve K(x) = b with GMRES, left preconditioned with the real factorization P
function preconditioned_gmres(K, b, x0, P)
    Pinv(v) = (w = P \ [real(v) imag(v)]; complex.(w[:, 1], w[:, 2]))
    x, info = KrylovKit.linsolve(
        v -> Pinv(K(v)),
        Pinv(b),
        x0,
        KrylovKit.GMRES(; krylovdim = 30, maxiter = 100, tol = 1e-10 * norm(x0)),
    )
    info.converged > 0 || @warn "GMRES did not converge" info.normres
    x
end
