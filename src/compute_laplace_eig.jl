"""
    compute_laplace_eig(mesh, pde, eiglim, neig_max)

Compute the Laplace eigenvalues, eigenfunctions and first order moments of products of pairs of eigenfunctions.
"""
function compute_laplace_eig(model, eiglim = Inf, neig_max = Inf)

    # Measure function evaluation time
    starttime = Base.time()

    # Extract parameters
    @unpack mesh, D, κ, T₂ = model
    ncompartment = length(mesh.points)

    # Assemble finite element matrices compartment-wise
    M_cmpts = []
    S_cmpts = []
    Mx_cmpts = [[] for dim = 1:3]
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt]
        facets = mesh.facets[icmpt, :]
        elements = mesh.elements[icmpt]
        volumes, _ = get_mesh_volumes(points, elements)

        # Assemble mass, stiffness and flux matrices
        push!(M_cmpts, assemble_mass_matrix(elements', volumes))
        push!(S_cmpts, assemble_stiffness_matrix(elements', points', D[icmpt]))

        # Assemble first order product moment matrices
        for dim = 1:3
            push!(Mx_cmpts[dim], assemble_mass_matrix(elements', volumes, points[dim, :]))
        end
    end

    # Assemble global finite element matrices
    M = blockdiag(M_cmpts...)
    S = blockdiag(S_cmpts...)
    R = blockdiag((M_cmpts ./ T₂)...)
    Mx = [blockdiag(Mx_cmpts[dim]...) for dim = 1:3]
    Q_blocks = assemble_flux_matrices(mesh.points, mesh.facets)
    Q = couple_flux_matrix(model, Q_blocks, false)

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
