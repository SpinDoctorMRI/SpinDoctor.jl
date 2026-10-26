"""
    solve(problem::Laplace)

Compute the Laplace eigenvalues, eigenfunctions and first order moments of products of pairs of eigenfunctions.

See [Li2020](@cite), [Agdestein2021](@cite).
"""
function solve(laplace::Laplace{T,dim}) where {T,dim}
    (; matrices, neig_max) = laplace
    (; M, S, R, Mx, Q) = matrices

    # Compute at most all eigenvalues in the given domain
    neig = Int(min(neig_max, size(M, 1)))

    @info join(
        [
            "Solving Laplace eigenvalue problem, computing $neig eigenvalues.",
            "Problem size: $(size(M, 1)) points.",
        ],
        "\n",
    )

    λ, ϕ = smallest_eigenpairs(S + Q, M, neig)

    # All Laplace eigenvalues are nonnegative (up to round-off)
    all(≥(-√eps(real(eltype(λ))) * maximum(abs, λ)), real.(λ)) ||
        @warn "Obtained negative eigenvalues for Laplace operator." findall(λ .< 0) λ[λ.<0]

    # Normalize eigenfunctions with mass weighting
    ϕ ./= sqrt.(sum(ϕ .* (M * ϕ), dims = 1))

    # Compute first order moments of product of pairs of eigenfunctions
    moments = [ϕ' * Mx[d] * ϕ for d = 1:dim]

    # Compute Laplace relaxation matrix
    massrelax = ϕ' * R * ϕ

    values = λ
    funcs = ϕ

    (; values, funcs, moments, massrelax)
end

"""
    λ, ϕ = smallest_eigenpairs(A, M, n)

Compute the `n` smallest eigenvalues `λ` and the corresponding `M`-orthonormal
eigenvectors `ϕ` of the generalized eigenproblem `A ϕ = λ M ϕ`, where `M` is symmetric
positive definite and `A` positive semi-definite.

The eigenvalues closest to a small negative shift `σ` (below the spectrum) are found with
Krylov iterations on the shift-and-invert operator `(A - σM)⁻¹ M`, which is self-adjoint in
the `M`-inner product when `A` is symmetric.
"""
function smallest_eigenpairs(A, M, n)
    N = size(M, 1)
    if N ≤ 1000 || n ≥ N ÷ 4
        # Small problem: dense decomposition
        Ad, Md = Matrix(A), Matrix(M)
        λ, ϕ = issymmetric(Ad) ? eigen(Symmetric(Ad), Symmetric(Md)) : eigen(Ad, Md)
        p = sortperm(real.(λ))[1:n]
        return λ[p], ϕ[:, p]
    end
    hermitian = issymmetric(A)
    σ = -1e-3 * tr(A) / tr(M)
    F = hermitian ? cholesky(Symmetric(A - σ * M)) : lu(A - σ * M)
    inner(x, y) = dot(x, M * y)
    op(v) = InnerProductVec(F \ (M * v.vec), inner)
    x₀ = InnerProductVec(1 .+ sin.(1:N), inner)
    ν, vecs, info = KrylovKit.eigsolve(
        op,
        x₀,
        n,
        :LR;
        ishermitian = hermitian,
        krylovdim = max(30, 2n + 10),
        tol = 1e-12,
    )
    info.converged ≥ n || @warn "Only $(info.converged) of $n eigenpairs converged"
    hermitian || (ν = real.(ν))
    λ = σ .+ 1 ./ ν[1:n]
    ϕ = reduce(hcat, [real.(v.vec) for v ∈ vecs[1:n]])
    p = sortperm(λ)
    λ[p], ϕ[:, p]
end
