"""
    eig2length(λ, D)

Convert Laplace eigenvalue to length scale.
"""
eig2length(λ, D) = λ > 0 ? π * √(D / λ) : Inf

"""
    length2eig(length, D)

Convert length scale to Laplace eigenvalue.
"""
length2eig(length, D) = D * (π / length)^2

"""
    limit_lengthscale(lap_eig, λ_max)

Only keep modes with length scales larger than minimum
"""
function limit_lengthscale(lap_eig, λ_max)
    (; values, funcs, moments, massrelax) = lap_eig
    inds = values .≤ λ_max
    values = values[inds]
    funcs = funcs[:, inds]
    moments = [m[inds, inds] for m ∈ moments]
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
        @warn "Obtained negative eigenvalues for Laplace operator." findall(λ .< 0) λ[λ .< 0]

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

"""
    solve_mf(model, matrices, lap_eig, gradient; ninterval = 500)

Compute the magnetization at the echo time using the matrix formalism, in the basis of the
Laplace eigenfunctions `lap_eig` (see [`solve_laplace`](@ref)). Gradients that are not
interval-wise constant are approximated by `ninterval` constant intervals.

See [Li2020](@cite), [Agdestein2021](@cite).
"""
function solve_mf(
    model::Model{TT,dim},
    matrices,
    lap_eig,
    gradient;
    ninterval = 500,
) where {TT,dim}
    (; γ) = model
    (; M) = matrices

    # Laplace eigenmodes
    λ = lap_eig.values
    ϕ = lap_eig.funcs
    Ax = lap_eig.moments
    T = lap_eig.massrelax
    neig = length(λ)

    ρ = initial_conditions(model)

    # Project initial spin density onto Laplace eigenfunction basis
    ν = ϕ' * (M * ρ)

    # Laplace operator in Laplace eigenfunction basis
    L = diagm(λ)

    # Bloch-Torrey operator in Laplace eigenfunction basis for given gradient
    function K!(K, g)
        if dim == 2
            @. K = L + T + im * γ * (g[1] * Ax[1] + g[2] * Ax[2])
        elseif dim == 3
            @. K = L + T + im * γ * (g[1] * Ax[1] + g[2] * Ax[2] + g[3] * Ax[3])
        end
    end

    K = zeros(eltype(ρ), neig, neig)

    if isa(gradient, ScalarGradient) && isconstant(gradient.profile)
        # The gradient is constant on each interval
        t = intervals(gradient.profile)
        for i = 1:(length(t)-1)
            δᵢ = t[i+1] - t[i]
            tᵢ = (t[i+1] + t[i]) / 2
            K!(K, gradient(tᵢ))
            expmv!(-δᵢ, K, ν)
        end
    else
        t = LinRange(0, echotime(gradient), ninterval + 1)
        for i = 1:ninterval
            δᵢ = t[i+1] - t[i]
            gᵢ = (gradient(t[i+1]) + gradient(t[i])) / 2
            K!(K, gᵢ)
            expmv!(-δᵢ, K, ν)
        end
    end

    # Final magnetization coefficients in finite element nodal basis
    ϕ * ν
end

"""
    compute_mf_diffusion_tensor(mesh, M, lap_eig, gradient::ScalarGradient)

Compute effective diffusion tensors using the matrix formalism.
"""
function compute_mf_diffusion_tensor(mesh, M, lap_eig, gradient::ScalarGradient)
    volumes = get_cmpt_volumes(mesh)
    p = reduce(hcat, mesh.points)
    a = p * M * lap_eig.funcs
    j = j_integral.(gradient.profile, lap_eig.values)
    D = a * Diagonal(j) * a' / sum(volumes)
end
