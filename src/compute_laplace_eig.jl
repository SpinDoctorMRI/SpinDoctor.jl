"""
    compute_laplace_eig(mesh, domain, eiglim, neig_max)

Compute the Laplace eigenvalues, eigenfunctions and first order moments of products of pairs of eigenfunctions.
"""
function compute_laplace_eig(mesh, domain, eiglim=Inf, neig_max=Inf)

    # Extract parameters
    @unpack ncompartment = mesh
    @unpack σ, κ, T₂ = domain

    # Assemble finite element matrices compartment-wise
    fem_mat_cmpts = (
        M = [],
        S = [],
        Q = [],
        Mx = [[] for idim=1:3]
    )
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt];
        facets = mesh.facets[icmpt, :];
        elements = mesh.elements[icmpt];
        volumes = get_mesh_volumes(points, elements);

        # Assemble mass, stiffness and flux matrices
        push!(fem_mat_cmpts.M, assemble_mass_matrix(elements', volumes))
        push!(fem_mat_cmpts.S, assemble_stiffness_matrix(elements', points', σ[icmpt]))
        push!(fem_mat_cmpts.Q, assemble_flux_matrix_cmpt(points, facets, κ))

        # Assemble first order product moment matrices
        for idim = 1:3
            push!(fem_mat_cmpts.Mx[idim], assemble_mass_matrix(elements', volumes, points[idim, :]))
        end
    end

    # Assemble global finite element matrices
    M = blockdiag(fem_mat_cmpts.M...)
    S = blockdiag(fem_mat_cmpts.S...)
    Q = couple_flux_matrix(mesh, fem_mat_cmpts.Q)
    Mx = [blockdiag(fem_mat_cmpts.Mx[idim]...) for idim = 1:3]

    # Compute at most all eigenvalues in the given domain
    neig = Int(min(neig_max, size(M, 1)));

    display("Solving Laplace eigenvalue problem, computing $neig eigenvalues.")
    display("Problem size: $(size(M, 1)) points.")

    # Solve generalized eigenproblem, computing the smallest eigenvalues only.
    # If 2 * neig_max_domain >= nnode, a full decomposition is performed,
    # calling the eig function inside eigs
    λ, ϕ = eigs(S+Q, M, nev=neig, which=:SR);
    # λ, ϕ = eigs(Hermitian(S+Q), Hermitian(M), nev=neig, which=:SR);
    # λ, ϕ = eigen(Hermitian(Matrix(S+Q)), Hermitian(Matrix(M)))

    # All Laplace eigenvalues are nonnegative
    all(0 .≤ λ) || @warn "Obtained negative eigenvalues for Laplace operator." findall(λ .< 0) λ[λ .< 0]

    # Only keep modes with length scales larger than minimum
    inds = λ .≤ eiglim
    λ = λ[inds]
    ϕ = ϕ[:, inds]
    isinf(eiglim) || length(λ) < neig || @warn "No eigenvalues were outside the interval. Consider increasing `neig_max`." eiglim neig_max

    # Normalize eigenfunctions with mass weighting
    ϕ ./= .√sum(ϕ .* (M * ϕ), dims=1)

    # Compute first order moments of product of pairs of eigenfunctions
    moments = cat([ϕ' * Mx[dim] * ϕ for dim = 1:3]..., dims=3)

    values = λ
    funcs = ϕ

    (; values, funcs, moments)
end
