"""
    solve(problem::Laplace)

Compute the Laplace eigenvalues, eigenfunctions and first order moments of products of pairs of eigenfunctions.
"""
function solve(laplace::Laplace)
    (; model, matrices, neig_max) = laplace
    (; mesh) = model
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

    # Solve generalized eigenproblem, computing the smallest eigenvalues only.
    # If 2 * neig_max_domain >= nnode, a full decomposition is performed,
    # calling the eig function inside eigs
    λ, ϕ = eigs(S + Q, M, nev = neig, which = :SR)
    # λ, ϕ = eigs(Hermitian(S + Q), Hermitian(M), nev = neig, which = :SR);
    # λ, ϕ = eigen(Hermitian(Matrix(S + Q)), Hermitian(Matrix(M)))

    # TODO: Consider KrylovKit.jl 

    # All Laplace eigenvalues are nonnegative
    all(0 .≤ λ) ||
        @warn "Obtained negative eigenvalues for Laplace operator." findall(λ .< 0) λ[λ.<0]

    # Normalize eigenfunctions with mass weighting
    ϕ ./= .√sum(ϕ .* (M * ϕ), dims = 1)

    # Compute first order moments of product of pairs of eigenfunctions
    moments = [ϕ' * Mx[dim] * ϕ for dim = 1:3]

    # Compute Laplace relaxation matrix
    massrelax = ϕ' * R * ϕ

    values = λ
    funcs = ϕ

    (; values, funcs, moments, massrelax)
end
