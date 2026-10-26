"""
    solve_analytical_mf(laplace, volumes, gradient)

Compute the signal in multilayered plates, cylinders or spheres using an analytical matrix
formalism solution, from the eigendecomposition `laplace` of
[`solve_analytical_laplace`](@ref) and the compartment `volumes`. Only implemented for
[`PGSE`](@ref) sequences.

This function is based on the following articles and corresponding code:
    [1] [Grebenkov2007](@cite)
    [2] [Grebenkov2010a](@cite)
"""
function solve_analytical_mf(laplace, volumes, gradient::ScalarGradient)
    (; r, ρ, γ, Λ, Br, B, U) = laplace

    gradient.profile isa PGSE || error("Only implemeted for PGSE")

    w = γ * gradient.amplitude * 1e12 * r[end]
    δ = gradient.profile.δ * 1e-6
    Δ = gradient.profile.Δ * 1e-6

    U₀ = complex(U)
    U = complex(U)

    expmv!(-δ, Λ + Br + im * w * B, U)
    expmv!(-(Δ - δ), Λ + Br, U)
    expmv!(-δ, Λ + Br - im * w * B, U)

    U₀'U * volumes'ρ
end
