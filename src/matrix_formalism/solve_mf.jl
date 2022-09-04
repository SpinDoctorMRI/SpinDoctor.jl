"""
    solve(problem::MatrixFormalism, gradient; ninterval = 500)

Solve for magnetization using Matrix Formalism.
"""
function solve(problem::MatrixFormalism{TT,dim}, gradient; ninterval = 500) where {TT,dim}
    (; model, matrices, lap_eig) = problem
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
        for i = 1:length(t)-1
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
