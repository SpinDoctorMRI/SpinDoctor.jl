# Integrals of time profiles for the short-time approximation and the matrix formalism

"""
    fresnels(x)

Compute the Fresnel S integral of type 0 (integral of `sin(pi t^2 / 2)`).
"""
fresnels(x) =
    (1 + im) / 4 * (erf((1 + im)sqrt(pi) / 2 * x) - im * erf((1 - im)sqrt(pi) / 2 * x))

"""
    fresnelc(x)

Compute the Fresnel C integral of type 0 (integral of `cos(pi t^2 / 2)`).
"""
fresnelc(x) =
    (1 - im) / 4 * (erf((1 + im)sqrt(pi) / 2 * x) + im * erf((1 - im)sqrt(pi) / 2 * x))

"""
    diffusion_time_sta(f)

Diffusion time of the time profile `f` in the short-time approximation of the ADC (see
[`compute_adc_sta`](@ref)).
"""
function diffusion_time_sta end

function diffusion_time_sta(f::PGSE)
    (; delta, Delta) = f
    out =
        4 / 35 * (
            (Delta + delta)^(7 / 2) + (Delta - delta)^(7 / 2) - 2Delta^(7 / 2) -
            2delta^(7 / 2)
        ) / (delta^2 * (Delta - delta / 3))
    out^2
end

function diffusion_time_sta(f::DoublePGSE)
    (; delta, Delta, tpause) = f
    tm = delta + tpause
    out =
        2 / 35 * (
            (2Delta + tm + delta)^(7 / 2) +
            (2Delta + tm - delta)^(7 / 2) +
            (tm + delta)^(7 / 2) +
            (tm - delta)^(7 / 2) - 2(2Delta + tm)^(7 / 2) - 2(Delta + tm + delta)^(7 / 2) -
            2(Delta + tm - delta)^(7 / 2) +
            2(Delta + delta)^(7 / 2) +
            2(Delta - delta)^(7 / 2) - 2tm^(7 / 2) + 4(Delta + tm)^(7 / 2) -
            4Delta^(7 / 2) - 4delta^(7 / 2)
        ) / (delta^2 * (Delta - delta / 3))
    out^2
end

function diffusion_time_sta(f::CosOGSE)
    (; delta, Delta) = f
    n = f.nperiod
    S = sin(2pi * n * Delta / delta)
    C = cos(2pi * n * Delta / delta)
    FS1 = fresnels(2sqrt(n))
    FS2 = fresnels(2sqrt(Delta * n / delta))
    FS3 = fresnels(2sqrt((Delta + delta) * n / delta))
    FS4 = fresnels(2sqrt((Delta - delta) * n / delta))
    FC1 = fresnelc(2sqrt(n))
    FC2 = fresnelc(2sqrt(Delta * n / delta))
    FC3 = fresnelc(2sqrt((Delta + delta) * n / delta))
    FC4 = fresnelc(2sqrt((Delta - delta) * n / delta))
    out =
        3 / 8 / sqrt(n * delta) * (
            2Delta * FC2 * C + delta * FC4 * C - delta * FC3 * C - Delta * FC3 * C -
            Delta * FC4 * C +
            2Delta * FS2 * S +
            delta * FS4 * S - delta * FS3 * S - Delta * FS3 * S - Delta * FS4 * S +
            2delta * FC1
        ) +
        9sqrt(delta) / 32pi / n^(3 / 2) *
        (+2FS2 * C - FS3 * C - FS4 * C - 2FC2 * S + FC3 * S + FC4 * S + 2FS1)
    out^2
end

function diffusion_time_sta(f::SinOGSE)
    (; delta, Delta) = f
    n = f.nperiod
    S = sin(2pi * n * Delta / delta)
    C = cos(2pi * n * Delta / delta)
    FS1 = fresnels(2sqrt(n))
    FS2 = fresnels(2sqrt(Delta * n / delta))
    FS3 = fresnels(2sqrt((Delta + delta) * n / delta))
    FS4 = fresnels(2sqrt((Delta - delta) * n / delta))
    FC1 = fresnelc(2sqrt(n))
    FC2 = fresnelc(2sqrt(Delta * n / delta))
    FC3 = fresnelc(2sqrt((Delta + delta) * n / delta))
    FC4 = fresnelc(2sqrt((Delta - delta) * n / delta))
    out =
        (
            64pi * n^(3 / 2) * Delta^(3 / 2) +
            64pi * n^(3 / 2) * delta^(3 / 2) +
            42delta^(3 / 2) * FS2 * C +
            42delta^(3 / 2) * FS1 - 42delta^(3 / 2) * FC2 * S +
            32pi * n^(3 / 2) * sqrt(Delta - delta) * delta -
            32pi * n^(3 / 2) * sqrt(Delta + delta) * Delta -
            32pi * n^(3 / 2) * sqrt(Delta + delta) * delta -
            32pi * n^(3 / 2) * sqrt(Delta - delta) * Delta +
            24pi * n * Delta * sqrt(delta) * FC2 * C +
            24pi * n * Delta * sqrt(delta) * FS2 * S +
            24pi * n * delta^(3 / 2) * FC1 +
            21delta^(3 / 2) * FC3 * S +
            21delta^(3 / 2) * FC4 * S - 21delta^(3 / 2) * FS3 * C -
            21delta^(3 / 2) * FS4 * C +
            12pi * n * delta^(3 / 2) * FS4 * S +
            12pi * n * delta^(3 / 2) * FC4 * C - 12pi * n * delta^(3 / 2) * FC3 * C -
            12pi * n * delta^(3 / 2) * FS3 * S - 12pi * n * Delta * sqrt(delta) * FC4 * C -
            12pi * n * Delta * sqrt(delta) * FC3 * C -
            12pi * n * Delta * sqrt(delta) * FS3 * S -
            12pi * n * Delta * sqrt(delta) * FS4 * S
        ) / (96pi * n^(3 / 2) * delta)
    out^2
end

"""
    j_integral(f, lambda)

Compute the integral

```math
j(f, \\lambda) = \\lambda \\frac{\\int_0^{T_\\text{echo}} F(t) \\int_0^t
\\mathrm{e}^{-\\lambda(t - s)} f(s) \\, \\mathrm{d} s \\, \\mathrm{d}
t}{\\int_0^{T_\\text{echo}} F^2(t) \\, \\mathrm{d} t}
```

for the time profile `f` and Laplace eigenvalue `lambda`.
"""
function j_integral end

function j_integral(f::PGSE, lambda)
    (; delta, Delta) = f

    if lambda < 1e-7
        # Use Taylor expansion when lambda is close to 0
        # to improve numerical stability
        lambda - lambda^2 * Delta^2 / 2(Delta - delta / 3) +
        lambda^3 * (10Delta^3 + 5Delta * delta^2 - delta^3) / 20(3 * Delta - delta) -
        lambda^4 * (Delta^4 + Delta^2 * delta^2) / 8(3Delta - delta) +
        lambda^5 * (21Delta^5 + 35Delta^3 * delta^2 + 7Delta * delta^4 - delta^5) /
        840(3Delta - delta)
    else
        -(
            exp(-lambda * (Delta + delta)) + exp(-lambda * (Delta - delta)) -
            2exp(-lambda * delta) - 2exp(-lambda * Delta) + 2(1 - lambda * delta)
        ) / (lambda^2 * delta^2 * (Delta - delta / 3))
    end
end

function j_integral(f::DoublePGSE, lambda)
    (; delta, Delta, tpause) = f
    tm = tpause + delta

    if lambda < 1e-7
        # Use Taylor expansion when lambda is close to 0
        # to improve numerical stability
        lambda - lambda^2 * 3Delta^2 / (3Delta - delta) +
        lambda^3 * (40Delta^3 + 5Delta * delta^2 - delta^3 + 30Delta^2 * tm) /
        20(3Delta - delta) -
        lambda^4 * (4Delta^4 + Delta^2 * delta^2 + 6Delta^3 * tm + 3Delta^2 * tm^2) /
        4(3Delta - delta) +
        lambda^5 * (
            336Delta^5 + 140Delta^3 * delta^2 + 7Delta * delta^4 - delta^5 +
            735Delta^4 * tm +
            105Delta^2 * delta^2 * tm +
            630Delta^3 * tm^2 +
            210Delta^2 * tm^3
        ) / 840(3Delta - delta)
    else
        -(
            exp(-lambda * (2Delta + tm - delta)) +
            exp(-lambda * (2Delta + tm + delta)) +
            exp(-lambda * (tm - delta)) +
            exp(-lambda * (tm + delta)) - 2exp(-lambda * (Delta + tm - delta)) -
            2exp(-lambda * (Delta + tm + delta)) - 2exp(-lambda * (2Delta + tm)) +
            2exp(-lambda * (Delta - delta)) +
            2exp(-lambda * (Delta + delta)) - 2exp(-lambda * tm) +
            4 +
            4exp(-lambda * (Delta + tm)) - 4exp(-lambda * Delta) - 4exp(-lambda * delta) -
            4lambda * delta
        ) / (2lambda^2 * delta^2 * (Delta - delta / 3))
    end
end

function j_integral(f::CosOGSE, lambda)
    (; delta, Delta, nperiod) = f
    n = nperiod

    if lambda < 1e-7
        # Use Taylor expansion when lambda is close to 0
        # to improve numerical stability
        lambda - lambda^3 * delta^2 * 3 / 4(n * pi)^2 +
        lambda^5 * (12Delta * n^2 * pi^2 * delta^3 + 15delta^4 - 4n^2 * pi^2 * delta^4) /
        (48n^4 * pi^4) - lambda^6 * Delta^2 * delta^3 / 8(n * pi)^2
    else
        4n^2 *
        pi^2 *
        lambda *
        (
            -exp(-lambda * (Delta + delta)) * lambda * delta -
            exp(-lambda * (Delta - delta)) * lambda * delta +
            2exp(-lambda * Delta) * lambda * delta +
            2exp(-lambda * delta) * lambda * delta +
            4n^2 * pi^2 +
            (lambda * delta - 2) * lambda * delta
        ) / (4n^2 * pi^2 + lambda^2 * delta^2)^2
    end
end

function j_integral(f::SinOGSE, lambda)
    (; delta, Delta, nperiod) = f
    n = nperiod

    if lambda < 1e-7
        # Use Taylor expansion when lambda is close to 0
        # to improve numerical stability
        lambda +
        lambda^3 * delta * (4n^2 * pi^2 * delta - 12Delta * n^2 * pi^2 - 15delta) /
        (36n^2 * pi^2) +
        lambda^4 * delta * Delta^2 / 6 +
        lambda^5 *
        delta *
        (
            120Delta * n^2 * pi^2 * delta^2 + 4n^4 * pi^4 * delta^3 -
            20Delta * n^4 * pi^4 * delta^2 + 105delta^3 - 40n^2 * pi^2 * delta^3 -
            40Delta^3 * n^4 * pi^4
        ) / (720n^4 * pi^4) +
        lambda^6 * (
            (Delta^4 * delta^2 + Delta^2 * delta^4) / 24 -
            (Delta^2 * delta^4) / (4n^2 * pi^2)
        ) / 3delta
    else
        4pi^2 *
        n^2 *
        (
            (lambda * delta)^3 +
            4n^2 *
            pi^2 *
            (
                exp(-lambda * (Delta + delta)) + exp(-lambda * (Delta - delta)) -
                2exp(-lambda * delta) - 2exp(-lambda * Delta) +
                lambda * delta +
                2
            )
        ) / (3delta * (4n^2 * pi^2 + lambda^2 * delta^2)^2)
    end
end
