# Matrix formalism: Laplace eigendecomposition and magnetization in the eigenfunction basis

"""
    eig2length(lambda, D)

Convert a Laplace eigenvalue to a length scale, for the diffusivity `D`.
"""
eig2length(lambda, D) = lambda > 0 ? pi * sqrt(D / lambda) : Inf

"""
    length2eig(length, D)

Convert a length scale to a Laplace eigenvalue, for the diffusivity `D`.
"""
length2eig(length, D) = D * (pi / length)^2

"""
    limit_lengthscale(lap_eig, lambda_max)

Keep the Laplace eigenmodes with eigenvalues up to `lambda_max` (length scales of at least
`eig2length(lambda_max, D)`).
"""
function limit_lengthscale(lap_eig, lambda_max)
    (; values, funcs, moments, massrelax) = lap_eig
    inds = values .<= lambda_max
    values = values[inds]
    funcs = funcs[:, inds]
    moments = [m[inds, inds] for m in moments]
    massrelax = massrelax[inds, inds]
    (; values, funcs, moments, massrelax)
end

"""
    solve_laplace(model, matrices; neig_max)

Compute the `neig_max` smallest Laplace eigenvalues `values`, the eigenfunctions `funcs`,
the first order moments of products of pairs of eigenfunctions `moments` and the
relaxation matrix `massrelax` in the eigenfunction basis.

See [Li2020](@cite), [Agdestein2021](@cite).
"""
function solve_laplace(model::Model{T,dim}, matrices; neig_max) where {T,dim}
    (; M, S, R, Mx, Q) = matrices

    # Compute at most all eigenvalues in the given domain
    neig = min(neig_max, size(M, 1))
    @info "Solving the Laplace eigenvalue problem for $neig eigenvalues ($(size(M, 1)) points)"

    lambda, phi = smallest_eigenpairs(S + Q, M, neig)

    # All Laplace eigenvalues are nonnegative (up to round-off)
    all(>=(-sqrt(eps(real(eltype(lambda)))) * maximum(abs, lambda)), real.(lambda)) ||
        @warn "Obtained negative eigenvalues for Laplace operator." findall(lambda .< 0) lambda[lambda .< 0]

    # Normalize eigenfunctions with mass weighting
    phi ./= sqrt.(sum(phi .* (M * phi), dims = 1))

    # Compute first order moments of product of pairs of eigenfunctions
    moments = [phi' * Mx[d] * phi for d = 1:dim]

    # Compute Laplace relaxation matrix
    massrelax = phi' * R * phi

    (; values = lambda, funcs = phi, moments, massrelax)
end

"""
    lambda, phi = smallest_eigenpairs(A, M, n)

Compute the `n` smallest eigenvalues `lambda` and the corresponding `M`-orthonormal
eigenvectors `phi` of the generalized eigenproblem `A phi = lambda M phi`, where `M` is symmetric
positive definite and `A` positive semi-definite.

The eigenvalues closest to a small negative shift `sigma` (below the spectrum) are found with
Krylov iterations on the shift-and-invert operator `(A - sigma M)^(-1) M`, which is
self-adjoint in the `M`-inner product when `A` is symmetric.
"""
function smallest_eigenpairs(A, M, n)
    N = size(M, 1)
    if N <= 1000 || n >= div(N, 4)
        # Small problem: dense decomposition
        Ad, Md = Matrix(A), Matrix(M)
        lambda, phi = issymmetric(Ad) ? eigen(Symmetric(Ad), Symmetric(Md)) : eigen(Ad, Md)
        p = sortperm(real.(lambda))[1:n]
        return lambda[p], phi[:, p]
    end
    hermitian = issymmetric(A)
    sigma = -1e-3 * tr(A) / tr(M)
    F = hermitian ? cholesky(Symmetric(A - sigma * M)) : lu(A - sigma * M)
    inner(x, y) = dot(x, M * y)
    op(v) = InnerProductVec(F \ (M * v.vec), inner)
    x0 = InnerProductVec(1 .+ sin.(1:N), inner)
    nu, vecs, info = KrylovKit.eigsolve(
        op,
        x0,
        n,
        :LR;
        ishermitian = hermitian,
        krylovdim = max(30, 2n + 10),
        tol = 1e-12,
    )
    info.converged >= n || @warn "Only $(info.converged) of $n eigenpairs converged"
    hermitian || (nu = real.(nu))
    lambda = sigma .+ 1 ./ nu[1:n]
    phi = reduce(hcat, [real.(v.vec) for v in vecs[1:n]])
    p = sortperm(lambda)
    lambda[p], phi[:, p]
end

"""
    solve_mf(model, matrices, lap_eig, gradient; ninterval = 500)

Compute the magnetization at the echo time using the matrix formalism, in the basis of the
Laplace eigenfunctions `lap_eig` (see [`solve_laplace`](@ref)). Gradients that are not
interval-wise constant are approximated by `ninterval` constant intervals.

See [Li2020](@cite), [Agdestein2021](@cite).
"""
function solve_mf(
    model::Model{T,dim},
    matrices,
    lap_eig,
    gradient;
    ninterval = 500,
) where {T,dim}
    (; gamma) = model
    (; M) = matrices

    # Laplace eigenmodes
    lambda = lap_eig.values
    phi = lap_eig.funcs
    Ax = lap_eig.moments
    relax = lap_eig.massrelax
    neig = length(lambda)

    # Project the initial spin density onto the Laplace eigenfunction basis
    xi0 = initial_conditions(model)
    nu = phi' * (M * xi0)

    # Laplace operator in Laplace eigenfunction basis
    L = diagm(lambda)

    # Bloch-Torrey operator in Laplace eigenfunction basis for given gradient
    function K!(K, g)
        if dim == 2
            @. K = L + relax + im * gamma * (g[1] * Ax[1] + g[2] * Ax[2])
        elseif dim == 3
            @. K = L + relax + im * gamma * (g[1] * Ax[1] + g[2] * Ax[2] + g[3] * Ax[3])
        end
    end

    K = zeros(eltype(xi0), neig, neig)

    if isa(gradient, ScalarGradient) && isconstant(gradient.profile)
        # The gradient is constant on each interval
        t = intervals(gradient.profile)
        for i = 1:(length(t)-1)
            dt = t[i+1] - t[i]
            t_i = (t[i+1] + t[i]) / 2
            K!(K, gradient(t_i))
            expmv!(-dt, K, nu)
        end
    else
        t = LinRange(0, echotime(gradient), ninterval + 1)
        for i = 1:ninterval
            dt = t[i+1] - t[i]
            g_i = (gradient(t[i+1]) + gradient(t[i])) / 2
            K!(K, g_i)
            expmv!(-dt, K, nu)
        end
    end

    # Final magnetization coefficients in finite element nodal basis
    phi * nu
end

"""
    compute_mf_diffusion_tensor(mesh, M, lap_eig, gradient::ScalarGradient)

Compute the effective diffusion tensor of the sample using the matrix formalism, from the
mass matrix `M` and the Laplace eigendecomposition `lap_eig`.
"""
function compute_mf_diffusion_tensor(mesh, M, lap_eig, gradient::ScalarGradient)
    volumes = compartment_volumes(mesh)
    p = reduce(hcat, mesh.points)
    a = p * M * lap_eig.funcs
    j = j_integral.(gradient.profile, lap_eig.values)
    a * Diagonal(j) * a' / sum(volumes)
end
