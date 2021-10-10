"""
    compute_laplace_eig(mesh, matrices, pde, eiglim, neig_max)

Compute the Laplace eigenvalues, eigenfunctions and first order moments of products of pairs of eigenfunctions.
"""
function compute_laplace_eig(model, matrices, eiglim = Inf, neig_max = Inf)

    # Measure function evaluation time
    starttime = Base.time()

    # Extract parameters
    @unpack mesh, D, T₂ = model
    @unpack M, S, R, Mx, Q = matrices
    ncompartment = length(mesh.points)

    # Compute at most all eigenvalues in the given domain
    neig = Int(min(neig_max, size(M, 1)))

    println("Solving Laplace eigenvalue problem, computing $neig eigenvalues.")
    println("Problem size: $(size(M, 1)) points.")

    # Solve generalized eigenproblem, computing the smallest eigenvalues only.
    # If 2 * neig_max_domain >= nnode, a full decomposition is performed,
    # calling the eig function inside eigs
    λ, ϕ = eigs(S + Q, M, nev = neig, which = :SR)
    # λ, ϕ = eigs(Hermitian(S + Q), Hermitian(M), nev = neig, which = :SR);
    # λ, ϕ = eigen(Hermitian(Matrix(S + Q)), Hermitian(Matrix(M)))

    # All Laplace eigenvalues are nonnegative
    all(0 .≤ λ) ||
        @warn "Obtained negative eigenvalues for Laplace operator." findall(λ .< 0) λ[λ.<0]

    # Only keep modes with length scales larger than minimum
    inds = λ .≤ eiglim
    λ = λ[inds]
    ϕ = ϕ[:, inds]
    isinf(eiglim) ||
        length(λ) < neig ||
        @warn "No eigenvalues were outside the interval. Consider increasing `neig_max`." eiglim neig_max

    # Normalize eigenfunctions with mass weighting
    ϕ ./= .√sum(ϕ .* (M * ϕ), dims = 1)

    # Compute first order moments of product of pairs of eigenfunctions
    moments = cat([ϕ' * Mx[dim] * ϕ for dim = 1:3]..., dims = 3)

    # Compute Laplace relaxation matrix
    massrelax = ϕ' * R * ϕ

    values = λ
    funcs = ϕ
    totaltime = Base.time() - starttime

    (; values, funcs, moments, massrelax, totaltime)
end
