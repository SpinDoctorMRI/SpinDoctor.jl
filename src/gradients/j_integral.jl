"""
    j_integral(f, λ)

Compute the integral

```math
j(f, \\lambda) = \\lambda \\frac{\\int_0^{T_\\text{echo}} F(t) \\int_0^t
\\mathrm{e}^{-\\lambda(t - s)} f(s) \\, \\mathrm{d} s \\, \\mathrm{d}
t}{\\int_0^{T_\\text{echo}} F^2(t) \\, \\mathrm{d} t}
```

for the time profile `f` and Laplace eigenvalue `λ`.
"""
function j_integral end

function j_integral(f::PGSE, λ)
    (; δ, Δ) = f

    if λ < 1e-7
        # Use Taylor expansion when λ is close to 0
        # to improve numerical stability
        λ - λ^2 * Δ^2 / 2(Δ - δ / 3) + λ^3 * (10Δ^3 + 5Δ * δ^2 - δ^3) / 20(3 * Δ - δ) -
        λ^4 * (Δ^4 + Δ^2 * δ^2) / 8(3Δ - δ) +
        λ^5 * (21Δ^5 + 35Δ^3 * δ^2 + 7Δ * δ^4 - δ^5) / 840(3Δ - δ)
    else
        -(
            exp(-λ * (Δ + δ)) + exp(-λ * (Δ - δ)) - 2exp(-λ * δ) - 2exp(-λ * Δ) +
            2(1 - λ * δ)
        ) / (λ^2 * δ^2 * (Δ - δ / 3))
    end
end

function j_integral(f::DoublePGSE, λ)
    #J Compute the quantity J(λ) for the sequence
    #   An analytical expression is available for the DoublePGSE sequence
    (; δ, Δ, tpause) = f
    tm = tpause + δ

    if λ < 1e-7
        # Use Taylor expansion when λ is close to 0
        # to improve numerical stability
        λ - λ^2 * 3Δ^2 / (3Δ - δ) +
        λ^3 * (40Δ^3 + 5Δ * δ^2 - δ^3 + 30Δ^2 * tm) / 20(3Δ - δ) -
        λ^4 * (4Δ^4 + Δ^2 * δ^2 + 6Δ^3 * tm + 3Δ^2 * tm^2) / 4(3Δ - δ) +
        λ^5 * (
            336Δ^5 + 140Δ^3 * δ^2 + 7Δ * δ^4 - δ^5 +
            735Δ^4 * tm +
            105Δ^2 * δ^2 * tm +
            630Δ^3 * tm^2 +
            210Δ^2 * tm^3
        ) / 840(3Δ - δ)
    else
        -(
            exp(-λ * (2Δ + tm - δ)) +
            exp(-λ * (2Δ + tm + δ)) +
            exp(-λ * (tm - δ)) +
            exp(-λ * (tm + δ)) - 2exp(-λ * (Δ + tm - δ)) - 2exp(-λ * (Δ + tm + δ)) -
            2exp(-λ * (2Δ + tm)) +
            2exp(-λ * (Δ - δ)) +
            2exp(-λ * (Δ + δ)) - 2exp(-λ * tm) +
            4 +
            4exp(-λ * (Δ + tm)) - 4exp(-λ * Δ) - 4exp(-λ * δ) - 4λ * δ
        ) / (2λ^2 * δ^2 * (Δ - δ / 3))
    end
end

function j_integral(f::CosOGSE, λ)
    #J Compute the quantity J(λ) for the sequence
    #   An analytical expression is available for the CosOGSE sequence
    (; δ, Δ, nperiod) = f
    n = nperiod

    if λ < 1e-7
        # Use Taylor expansion when λ is close to 0
        # to improve numerical stability
        λ - λ^3 * δ^2 * 3 / 4(n * π)^2 +
        λ^5 * (12Δ * n^2 * π^2 * δ^3 + 15δ^4 - 4n^2 * π^2 * δ^4) / (48n^4 * π^4) -
        λ^6 * Δ^2 * δ^3 / 8(n * π)^2
    else
        4n^2 *
        π^2 *
        λ *
        (
            -exp(-λ * (Δ + δ)) * λ * δ - exp(-λ * (Δ - δ)) * λ * δ +
            2exp(-λ * Δ) * λ * δ +
            2exp(-λ * δ) * λ * δ +
            4n^2 * π^2 +
            (λ * δ - 2) * λ * δ
        ) / (4n^2 * π^2 + λ^2 * δ^2)^2
    end
end

function j_integral(f::SinOGSE, λ)
    #J Compute the quantity J(λ) for the sequence
    #   An analytical expression is available for the SinOGSE sequence
    (; δ, Δ, nperiod) = f
    n = nperiod

    if λ < 1e-7
        # Use Taylor expansion when λ is close to 0
        # to improve numerical stability
        λ +
        λ^3 * δ * (4n^2 * π^2 * δ - 12Δ * n^2 * π^2 - 15δ) / (36n^2 * π^2) +
        λ^4 * δ * Δ^2 / 6 +
        λ^5 *
        δ *
        (
            120Δ * n^2 * π^2 * δ^2 + 4n^4 * π^4 * δ^3 - 20Δ * n^4 * π^4 * δ^2 + 105δ^3 -
            40n^2 * π^2 * δ^3 - 40Δ^3 * n^4 * π^4
        ) / (720n^4 * π^4) +
        λ^6 * ((Δ^4 * δ^2 + Δ^2 * δ^4) / 24 - (Δ^2 * δ^4) / (4n^2 * π^2)) / 3δ
    else
        4π^2 *
        n^2 *
        (
            (λ * δ)^3 +
            4n^2 *
            π^2 *
            (
                exp(-λ * (Δ + δ)) + exp(-λ * (Δ - δ)) - 2exp(-λ * δ) - 2exp(-λ * Δ) +
                λ * δ +
                2
            )
        ) / (3δ * (4n^2 * π^2 + λ^2 * δ^2)^2)
    end
end
