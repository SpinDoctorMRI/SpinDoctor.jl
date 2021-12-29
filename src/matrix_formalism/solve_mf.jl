"""
    solve(problem::MatrixFormalism, matrices, lap_eig, gradient)

Solve for magnetization using Matrix Formalism.
"""
function solve(problem::MatrixFormalism, gradient)
    (; model, matrices, lap_eig, ninterval) = problem
    (; mesh, ρ, γ) = model
    (; M, M_cmpts) = matrices

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
    function K!(K, g⃗)
        @. K = L + T + im * γ * (g⃗[1] * Ax[1] + g⃗[2] * Ax[2] + g⃗[3] * Ax[3])
    end
    
    K = zeros(eltype(ρ), neig, neig)

    if isa(gradient, ScalarGradient) && isconstant(gradient.profile)
        # The gradient is constant on each interval
        t = intervals(gradient.profile)
        for i = 1:length(t)-1
            δᵢ = t[i + 1] - t[i]
            tᵢ = (t[i + 1] + t[i]) / 2
            K!(K, gradient(tᵢ))
            expmv!(-δᵢ, K, ν)
        end
    else
        t = LinRange(0, echotime(gradient), ninterval + 1)
        for i = 1:ninterval
            δᵢ = t[i + 1] - t[i]
            g⃗ᵢ = (gradient(t[i + 1]) + gradient(t[i])) / 2
            K!(K, g⃗ᵢ)
            expmv!(-δᵢ, K, ν)
        end
    end

    # Final magnetization coefficients in finite element nodal basis
    ϕ * ν
end
