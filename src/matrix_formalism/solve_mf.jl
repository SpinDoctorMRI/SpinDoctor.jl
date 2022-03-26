"""
    solve(problem::MatrixFormalism, gradient; ninterval = 500)

Solve for magnetization using Matrix Formalism.
"""
function solve(problem::MatrixFormalism{T}, gradient; ninterval = 500) where{T}
    (; model, matrices, lap_eig) = problem
    (; mesh, γ, T₂, ncompartment) = model
    (; M_cmpts, M) = matrices
    (; λ, ϕ, Ax, MT₂,ν) = lap_eig

    if isempty(MT₂)
        expT₂ = exp(-echotime(gradient)/first(T₂))
    end

    ξ_cmpts = VectorOfArrays{Complex{T}, 1}()
    for icmpt = 1:ncompartment
        ξ_tmp = Vector{Complex{T}}(undef,size(ϕ[icmpt],1))
        if isempty(MT₂)
            ν_tmp = evolve_laplace_coef( ν[icmpt], γ, gradient, Diagonal(λ[icmpt]), Ax[icmpt], ninterval)
            # Final magnetization coefficients in finite element nodal basis
            mul!(ξ_tmp, ϕ[icmpt], ν_tmp )
            ξ_tmp .= expT₂ .* ξ_tmp
        else
            ν_tmp = evolve_laplace_coef( ν[icmpt], γ, gradient, Diagonal(λ[icmpt]) + MT₂, Ax[icmpt], ninterval)
            # Final magnetization coefficients in finite element nodal basis
            mul!(ξ_tmp, ϕ[icmpt], ν_tmp )
        end
        push!(ξ_cmpts, ξ_tmp)
    end
    ξ = flatview(ξ_cmpts)
    S = compute_signal(M, ξ)
    S_cmpts = compute_signal.(M_cmpts, split_field(mesh, ξ))

    S, S_cmpts, ξ
end

function evolve_laplace_coef( ν, γ, gradient, LMT₂, Ax, ninterval)
    # Bloch-Torrey operator in Laplace eigenfunction basis for given gradient 
    neig = length(ν)
    K = zeros(eltype(ν), neig, neig)

    function K!(K, g⃗)
        @. K = LMT₂ + im * γ * (g⃗[1] * Ax[1] + g⃗[2] * Ax[2] + g⃗[3] * Ax[3])
    end

    if isa(gradient.profile, PGSE)
        K!(K, gradient( 0 ))
        # expmv!(-gradient.profile.δ, K, ν)
        ν_tmp, = exponentiate(K,-gradient.profile.δ,ν)
        if isdiag(LMT₂)
            diagLMT₂ = diag(LMT₂)
            lmul!(-(gradient.profile.Δ - gradient.profile.δ),diagLMT₂)
            @. ν_tmp = ν_tmp * exp(diagLMT₂)
        else
            ν_tmp, = exponentiate(LMT₂,-(gradient.profile.Δ - gradient.profile.δ), ν_tmp)
        end
        conj!(K)
        # expmv!(-gradient.profile.δ, K, ν)
        ν_tmp, = exponentiate(K,-gradient.profile.δ,ν_tmp)

    else
        ν_tmp = copy(ν)
        if isa(gradient, ScalarGradient) && isconstant(gradient.profile)
            # The gradient is constant on each interval
            t = intervals(gradient.profile)
            for i = 1:length(t)-1
                δᵢ = t[i+1] - t[i]
                tᵢ = (t[i+1] + t[i]) / 2
                K!(K, gradient(tᵢ))
                if isdiag(K)
                    ν_tmp = exp.(-δᵢ*diag(K) ) .* ν_tmp
                else
                    # expmv!(-δᵢ, K, ν)
                    ν_tmp, = exponentiate(K,-δᵢ,ν_tmp)
                end
            end
        else
            t = LinRange(0, echotime(gradient), ninterval + 1)
            for i = 1:ninterval
                δᵢ = t[i+1] - t[i]
                g⃗ᵢ = (gradient(t[i+1]) + gradient(t[i])) / 2
                K!(K, g⃗ᵢ)
                # expmv!(-δᵢ, K, ν)
                ν_tmp, = exponentiate(K,-δᵢ,ν_tmp)
            end
        end
    end
    ν_tmp
end