"""
    solve(problem::AnalyticalMatrixFormalism, gradient)

Compute the signal in a multilayered cylinder or sphere using an analytical
matrix formalism solution.

This function is based on the following articles and corresponding code:
    [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
        Rev. Mod.Phys. 79, 1077 (2007)
    [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted diffusion in
        multilayered structures, J. Magn. Reson. 205, 181-195 (2010).
"""
function solve(problem::AnalyticalMatrixFormalism, gradient::ScalarGradient)
    (; analytical_laplace, lap_mat, volumes) = problem
    (; r, ρ, γ) = analytical_laplace
    (; Λ, Br, B, U) = lap_mat

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
