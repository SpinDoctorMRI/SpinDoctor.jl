"""
    solve(problem::Laplace)

Compute the Laplace eigenvalues, eigenfunctions and first order moments of products of pairs of eigenfunctions.
"""
function solve(laplace::Laplace{T};rerun::Bool = false,savepath::NamedTuple = NamedTuple()) where{T}
    (; model, matrices, neig_max, length_scale) = laplace
    (; κ, T₂, D_avg) = model

    λ_max = length2eig(length_scale, D_avg)

    ## try to load saved results if it exists.
    if !isempty(savepath) && !rerun
        (; path, checksum) = savepath
        if ispath(path)
            filename = joinpath(path, "lap_eig_values_lengthscale-$(floor(Int,length_scale*1e4))-neigmax-$(neig_max)-$(checksum).jld2")
            if isfile(filename)
                λ,ls,Ax,MT₂,ν = load(filename, "λ", "ls","Ax","MT₂","ν")
                @info "Load eig_lap from "*filename
                filename = joinpath(path, "lap_eig_func_lengthscale-$(floor(Int,length_scale*1e4))-neigmax-$(neig_max)-$(checksum).jld2")
                ϕ = load(filename, "ϕ")
                return (; λ, ϕ, ls, Ax, MT₂, ν)
            else
                searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

                fileslist = searchdir(path,Regex("lap_eig_values_lengthscale-.*(-neigmax-)+.*-$(checksum).jld2" ) )
                if !isempty(fileslist)
                    for filename ∈ fileslist
                        actual_ls = parse(Int64,split(filename,'-')[2])
                        neigmax = parse(Int64,split(filename,'-')[4])

                        if length_scale ≥ actual_ls/1e4 && neig_max ≤ neigmax
                            filename = joinpath(path, filename)
                            λ,ls,Ax,MT₂,ν = load(filename, "λ", "ls","Ax","MT₂","ν")
                            @info "Load eig_lap from "*filename
                            filename = replace(filename, "lap_eig_values_lengthscale" => "lap_eig_func_lengthscale-")
                            ϕ = load(filename, "ϕ")
                            # truncation
                            ncompartment = length(λ)
                            for icmpt =  1:ncompartment
                                inds = λ[icmpt] .≤ λ_max
                                λ[icmpt] = λ[icmpt][inds]
                                ϕ[icmpt] = ϕ[icmpt][:, inds]
                                ls[icmpt] = ls[icmpt][inds]
                                Ax[icmpt] = [m[inds, inds] for m ∈ Ax[icmpt]]
                                if !all(y -> y == first(T₂), T₂)
                                    MT₂[icmpt] = MT₂[icmpt][inds, inds]
                                end
                                ν[icmpt] = ν[icmpt][inds]
                            end
                            return (; λ, ϕ, ls, Ax, MT₂, ν)
                        end
                    end
                end
            end
        end
    end

    ps = MKLPardisoSolver()
    Ax = VectorOfArrays{Matrix{T}, 1}()
    ϕ = VectorOfArrays{T, 2}()
    λ = VectorOfArrays{T, 1}()
    ls = VectorOfArrays{T, 1}()
    ν = VectorOfArrays{Complex{T}, 1}()

    if iszero(κ)
        (; ncompartment) = model
        (; M_cmpts, S_cmpts, Mx_cmpts, ρ_cmpts) = matrices

        actual_ls_cmpts = Vector{T}(undef, ncompartment)
        for icmpt = 1:ncompartment
            # eigen values problem
            values,funcs = eigendecompose(S_cmpts[icmpt], M_cmpts[icmpt], laplace, λ_max, ps)
            push!(λ, values)
            push!(ϕ, funcs)
            # Compute length scale
            ls_cmpts = eig2length.(values, D_avg)
            push!(ls, ls_cmpts)
            # Compute first order moments of product of pairs of eigenfunctions
            push!(Ax, [funcs' * Mx_cmpts[dim][icmpt] * funcs for dim = 1:3])
            
            # Project initial spin density onto Laplace eigenfunction basis
            push!(ν, funcs' * (transpose(M_cmpts[icmpt]) * ρ_cmpts[icmpt]) )

            actual_ls_cmpts[icmpt] = length(values) < size(M_cmpts[icmpt],1) ? length_scale : minimum(ls_cmpts)
        end
        # get actual length scale
        actual_ls = maximum(actual_ls_cmpts)
        MT₂ = Matrix{T}[]
    else
        (; M, S, R, Q, Mx, ρ_init) = matrices;
        # decide use Lanczos or Arnoldi
        isQherm = ishermitian(Q)
        # eigen values problem
        values,funcs = eigendecompose(S+Q, M, laplace, λ_max, ps; isQherm=isQherm)
        push!(λ, values)
        push!(ϕ, funcs)
        # Compute length scale
        ls_cmpts = eig2length.(values, D_avg)
        push!(ls, ls_cmpts)
        # Compute first order moments of product of pairs of eigenfunctions
        push!(Ax, [funcs' * Mx[dim] * funcs for dim = 1:3])
        
        if all(y -> y == first(T₂), T₂)
            MT₂ = Matrix{T}[]
        else
            # Compute Laplace relaxation matrix
            MT₂ = funcs' * R * funcs
        end

        # Project initial spin density onto Laplace eigenfunction basis
        push!(ν, funcs' * (transpose(M) * ρ_init) )

        # get actual length scale
        actual_ls = length(values) < size(M,1) ? length_scale : minimum(ls_cmpts)
    end

    ## try to save results if savepath is provided
    if !isempty(savepath)
        (; path, checksum) = savepath
        ispath(path) || mkpath(path)
        filename = joinpath(path, "lap_eig_values_lengthscale-$(floor(Int,actual_ls*1e4))-neigmax-$(neig_max)-$(checksum).jld2")
        jldsave(filename; λ,ls,Ax,MT₂,ν)
        filename = joinpath(path, "lap_eig_func_lengthscale-$(floor(Int,actual_ls*1e4))-neigmax-$(neig_max)-$(checksum).jld2")
        jldsave(filename; ϕ)
        @info "saved lapeig"
    end

    (; λ, ϕ, ls, Ax, MT₂, ν)
end

function eigendecompose(SQ, M, laplace::Laplace, λ_max, ps; isQherm::Bool=true)
    (;neig_max,ncv,tol,maxiter) = laplace;

    sizeofM = size(M, 1)

    # Compute at most all eigenvalues in the given domain
    neig = Int(min(neig_max, sizeofM))

    @info join(
        [
            "Solving $(isQherm ? "symmetric" : "unsymmetric") Laplace eigenvalue problem, computing $neig eigenvalues.",
            "Problem size: $(sizeofM) points.",
        ],
        "\n",
    )

    # Solve generalized eigenproblem, computing the smallest eigenvalues only.
    # If 2 * neig_max_domain >= nnode, a full decomposition is performed
    if sizeofM < 4e3 || sizeofM < 2*neig  # for small size, consider use eigen

        F = @time eigen(Hermitian(Matrix(SQ)), Hermitian(Matrix(M)))
        λ, ϕ = F

        sort_λ = zeros(Int, length(λ))
        sortperm!(sort_λ, λ)
        ϕ = ϕ[:, sort_λ]
        λ = λ[sort_λ]

    else
        if isQherm
            # set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
            lm, lm_phi, ps = @time construct_linear_map(SQ+1e-16.*M, transpose(M), ps)

            λ, ϕ = @time eigs(lm, nev = neig, ncv = ncv, which = :LR, tol = tol, maxiter = maxiter)
            ϕ = @time reduce(hcat,[lm_phi(x) for x ∈ eachcol(ϕ) ]);

            # TODO: Consider KrylovKit.jl
            # λ, ϕ, = @time eigsolve(lm, ones(sizeofM), neig, :LR,
            #     Lanczos(krylovdim = ncv, maxiter = maxiter, tol = 1e-16, verbosity = 1) )
            # ϕ = @time reduce(hcat,[lm_phi(x) for x ∈ ϕ ])
             
            λ .= -1e-16 .+ (1. ./ λ)

        else
            lm, ps = @time construct_nonsymmetric_linear_map(SQ+1e-20.*M, transpose(M), ps)
            λ, ϕ = @time eigs(lm, nev = neig, ncv = ncv, which = :LR, tol = tol, maxiter = maxiter)

            # λ, ϕ, = @time eigsolve(lm, ones(sizeofM), neig, :LR,
            #     Arnoldi(krylovdim = ncv, maxiter = maxiter,tol = 1e-14, verbosity = 0) )
            # ϕ = @time reduce(hcat,ϕ)


            λ .= -1e-20 .+ (1. ./ λ)
        end
        set_phase!(ps, Pardiso.RELEASE_ALL);
        pardiso(ps);
    end

    if abs(λ[1]) ≤ eps()  λ[1] = 0.0 end 
    # All Laplace eigenvalues are nonnegative
    all(≥(0), λ) || @warn "Obtained negative eigenvalues for Laplace operator." findall(λ .< 0) λ[λ.<0] 
    # if any(<(0), λ)
    #     index_neg = findall(λ .< 0)
    #     @warn "Obtained negative eigenvalues for Laplace operator." index_neg λ[index_neg]
    #     λ[index_neg] .= 0.0
    # end

    # truncation
    if !isinf(λ_max)
        inds = λ .≤ λ_max
        if all(inds)
            if neig == sizeofM
                @info "All eigenvalues were inside the interval. No truncation."
            else
                @info "No eigenvalues were outside the interval. Consider increasing neig_max."
            end
        else
            λ = λ[inds]
            ϕ = ϕ[:, inds]
        end
    else
        @info "No truncation"
    end

    # Normalize eigenfunctions with mass weighting
    ϕ ./= .√sum(ϕ .* (M * ϕ), dims = 1)

    λ,ϕ
end

function construct_linear_map(H, S,ps) 
    if eltype(H)<:AbstractFloat
        set_matrixtype!(ps, Pardiso.REAL_SYM_POSDEF)
    else
        set_matrixtype!(ps, Pardiso.COMPLEX_HERM_POSDEF)
    end
    set_nprocs!(ps, Threads.nthreads() )
    pardisoinit(ps)
    fix_iparm!(ps, :N)

    H_pardiso = get_matrix(ps, H, :N)
    b = rand(Float64, size(H, 1))

    set_phase!(ps, Pardiso.ANALYSIS_NUM_FACT)
    pardiso(ps, H_pardiso, b)

    return (
        LinearMap{eltype(H)}(
            (y, x) -> begin
                set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE_ONLY_BACKWARD)
                pardiso(ps, y, H_pardiso, x)
                set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE_ONLY_FORWARD)
                pardiso(ps, y, H_pardiso, S * y)
            end,
            size(H, 1);
            issymmetric=true,
            ismutating=true
        ),
        LinearMap{eltype(H)}(
            (y, x) -> begin
                set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE_ONLY_BACKWARD)
                pardiso(ps, y, H_pardiso, x)
            end,
            size(H, 1);
            ismutating=true
        ),
        ps,
    )
end

function construct_nonsymmetric_linear_map(H, S,ps)
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    if eltype(H)<:AbstractFloat
        set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    else
        set_matrixtype!(ps, Pardiso.COMPLEX_NONSYM)
    end
    set_nprocs!(ps, Threads.nthreads() )
    pardisoinit(ps)
    fix_iparm!(ps, :N)

    H_pardiso = get_matrix(ps, H, :N)
    b = rand(Float64, size(H, 1))

    set_phase!(ps, Pardiso.ANALYSIS_NUM_FACT)
    pardiso(ps, H_pardiso, b)

    return (
        LinearMap{eltype(H)}(
            (y, x) -> begin
                set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
                pardiso(ps, y, H_pardiso, S * x)
            end,
            size(H, 1);
            issymmetric=false,
            ismutating=true
        ),
        ps,
    )
end