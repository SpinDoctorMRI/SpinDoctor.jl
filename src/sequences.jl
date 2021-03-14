
"""
General time profile type (gradient sequence).

A `TimeProfile` should implement a calling method returning the value of the
time profile at time `t`. It can optionally overwrite the `echotime` `integral`,
and `bvalue_no_q` methods, the first providing a default value and the latter
two computing numerical approximation to the integral quantities if not
overwritten.
"""
abstract type TimeProfile end


"""
    f = PGSE(δ, Δ)

Pulsed Gradient Spin Echo sequence with pulse duration `δ` and time between
pulses `Δ-δ`.
"""
struct PGSE <: TimeProfile
    δ::Float64
    Δ::Float64
end


"""
    f = CosOGSE(δ, Δ, nperiod)

Oscillating Gradient Spin Echo sequence with two cos-pulses of duration `δ`
separated by a pause of duration `Δ-δ` for `nperiod` periods per pulse.
"""
struct CosOGSE <: TimeProfile
    δ::Float64
    Δ::Float64
    nperiod::Int
end


"""
    f = SinOGSE(δ, Δ, nperiod)

Oscillating Gradient Spin Echo sequence with two sin-pulses of duration `δ`
separated by a pause of duration `Δ-δ` for `nperiod` periods per pulse."""
struct SinOGSE <: TimeProfile
    δ::Float64
    Δ::Float64
    nperiod::Int
end


"""
    f = DoublePGSE(δ, Δ)

Double Pulsed Gradient Spin Echo sequence with four pulses of duration `δ`,
separated by pauses of duration `Δ-δ`, `0`, and `Δ-δ` repsectively.
"""
struct DoublePGSE <: TimeProfile
    δ::Float64
    Δ::Float64
end


"""
    ft = (f::TimeProfile)(t)

Return the time profile value at time `t`.
"""
function (f::TimeProfile)(t)
    error("Not implemented")
end

function (f::PGSE)(t)
    (t < f.δ) - (f.Δ ≤ t)
end

function (f::CosOGSE)(t)
    (
        (t < f.δ) * cos(2π * f.nperiod * t / f.δ) -
        (f.Δ ≤ t) * cos(2π * f.nperiod * (t - f.Δ) / f.δ)
    )
end

function (f::SinOGSE)(t)
    (
        (t < f.δ) * sin(2π * f.nperiod * t / f.δ) -
        (f.Δ ≤ t) * sin(2π * f.nperiod * (t - f.Δ) / f.δ)
    )
end

function (f::DoublePGSE)(t)
    ((t < f.δ) - (f.Δ ≤ t < f.Δ + f.δ) + (f.Δ + f.δ ≤ t < 2f.Δ + f.δ) - (2f.Δ + f.δ ≤ t))
end


"""
    TE = echotime(f)

Get echo time `TE` of time profile `f`. By default, it returns `f.Δ+f.δ`.

For the double gradient spin echo (`DoublePGSE`) sequence, the echo time is
doubled.
"""
function echotime(f::TimeProfile)
    f.Δ + f.δ
end

function echotime(f::DoublePGSE)
    2f.Δ + 2f.δ
end


"""
    integral(f, t)
    integral(f) -> integral(f, echotime(f))

Integral of time profile `f` between `0` and `t`. Unless specified, the echotime
is used as the upper integral limit.

For the `PGSE`, `SinOGSE`, `CosOGSE` and `DoublePGSE` sequences, analytical
expressions are available. Otherwise a numerical integral is computed.
"""
function integral(f::TimeProfile, t = echotime(f))
    quadgk(f, 0, t)
end

function integral(f::PGSE, t = echotime(f))
    ((t < f.δ) * t + (f.δ ≤ t) * f.δ - (f.Δ ≤ t) * (t - f.Δ))
end

function integral(f::CosOGSE, t = echotime(f))
    δ, Δ, n = f.δ, f.Δ, f.nperiod
    ((t < δ) * sin(2π * n * t / δ) - (Δ ≤ t) * sin(2π * n * (t - Δ) / δ)) * δ / (2π * n)
end

function integral(f::SinOGSE, t = echotime(f))
    δ, Δ, n = f.δ, f.Δ, f.nperiod
    ((t < δ) * (1 - cos(2π * n * t / δ)) - (Δ ≤ t) * (1 - cos(2π * n * (t - Δ) / δ))) * δ /
    (2π * n)
end

function integral(f::DoublePGSE, t = echotime(f))
    δ, Δ = f.δ, f.Δ
    tmid = Δ + δ
    (
        (t < δ) * t + (δ ≤ t < tmid) * δ - (Δ ≤ t < tmid) * (t - Δ) +
        (tmid ≤ t < tmid + δ) * (t - tmid) +
        (tmid + δ ≤ t) * δ - (tmid + Δ ≤ t) * (t - (tmid + Δ))
    )
end


"""
    b = bvalue_no_q(f)

Compute the time profile contribution to the b-value. To obtain the b-value,
multiply the result by q^2 = (γg)^2, where γ is the gyromagnetic ratio of the
water proton, g is the gradient amplitude, and b(q, f) = q^2 * bvalue_no_q(f).

"""
function bvalue_no_q(f::TimeProfile)
    quadgk(s -> integral(f, s)^2, 0, echotime(f))
end

function bvalue_no_q(f::PGSE)
    f.δ^2 * (f.Δ - f.δ / 3)
end

function bvalue_no_q(f::CosOGSE)
    f.δ^3 / (2π * f.nperiod)^2
end

function bvalue_no_q(f::SinOGSE)
    3f.δ^3 / (2π * f.nperiod)^2
end

function bvalue_no_q(f::DoublePGSE)
    2 * f.δ^2 * (f.Δ - f.δ / 3)
end
