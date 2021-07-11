
"""
General time profile type (gradient sequence).

A `TimeProfile` should implement a calling method returning the value of the
time profile at time `t`. It can optionally overwrite the `echotime` `integral`,
and `bvalue_no_q` methods, the first providing a default value and the latter
two computing numerical approximation to the integral quantities if not
overwritten.
"""
abstract type TimeProfile{T<:Real} <: Function end

"""
    f = PGSE(δ, Δ)

Pulsed Gradient Spin Echo sequence with pulse duration `δ` and time between
pulses `Δ-δ`.
"""
struct PGSE{T} <: TimeProfile{T}
    δ::T
    Δ::T
    function PGSE(δ::T, Δ::T) where {T}
        @assert δ ≤ Δ
        new{T}(δ, Δ)
    end
end
PGSE(δ::Int, Δ::Int) = PGSE(float(δ), float(Δ))
PGSE(δ::Int, Δ) = PGSE(promote(δ, Δ)...)
PGSE(δ, Δ::Int) = PGSE(promote(δ, Δ)...)


"""
    f = CosOGSE(δ, Δ, nperiod)

Oscillating Gradient Spin Echo sequence with two cos-pulses of duration `δ`
separated by a pause of duration `Δ-δ` for `nperiod` periods per pulse.
"""
struct CosOGSE{T} <: TimeProfile{T}
    δ::T
    Δ::T
    nperiod::Int
    function CosOGSE(δ::T, Δ::T, n) where {T}
        @assert δ ≤ Δ
        @assert n > zero(n)
        new{T}(δ, Δ, n)
    end
end

"""
    f = SinOGSE(δ, Δ, nperiod)

Oscillating Gradient Spin Echo sequence with two sin-pulses of duration `δ`
separated by a pause of duration `Δ-δ` for `nperiod` periods per pulse.
"""
struct SinOGSE{T} <: TimeProfile{T}
    δ::T
    Δ::T
    nperiod::Int
    function SinOGSE(δ::T, Δ::T, n) where {T}
        @assert δ ≤ Δ
        @assert n > zero(n)
        new{T}(δ, Δ, n)
    end
end

"""
    f = DoublePGSE(δ, Δ)

Double Pulsed Gradient Spin Echo sequence with four pulses of duration `δ`,
separated by pauses of duration `Δ-δ`, `0`, and `Δ-δ` repsectively.
"""
struct DoublePGSE{T} <: TimeProfile{T}
    δ::T
    Δ::T
    tpause::T
    function DoublePGSE(δ::T, Δ::T, tpause::T) where {T}
        @assert δ ≤ Δ
        @assert tpause ≥ zero(tpause)
        new{T}(δ, Δ, tpause)
    end
end



"""
    ft = (f::TimeProfile)(t)

Return the time profile value at time `t`.
"""
function (f::TimeProfile)(t)
    error("Not implemented")
end

function (f::PGSE)(t)
    (zero(t) ≤ t < f.δ) - (f.Δ ≤ t < f.Δ + f.δ)
end

function (f::DoublePGSE)(t)
    δ, Δ, p = f.δ, f.Δ, f.tpause
    (
        (zero(t) ≤ t < δ) - (Δ ≤ t < Δ + δ) + (p + Δ + δ ≤ t < p + Δ + 2δ) -
        (p + 2Δ + δ ≤ t < p + 2Δ + 2δ)
    )
end

function (f::CosOGSE)(t)
    δ, Δ, n = f.δ, f.Δ, f.nperiod
    (t < δ) * cos(2π * nperiod * t / δ) - (Δ ≤ t) * cos(2π * nperiod * (t - Δ) / δ)
end

function (f::SinOGSE)(t)
    δ, Δ, n = f.δ, f.Δ, f.nperiod
    (t < δ) * sin(2π * nperiod * t / δ) - (Δ ≤ t) * sin(2π * nperiod * (t - Δ) / δ)
end


"""
    integral(f, t = echotime(f))

Integral of time profile `f` between `0` and `t`. Unless specified, the echotime
is used as the upper integral limit.

For the `PGSE`, `SinOGSE`, `CosOGSE` and `DoublePGSE` sequences, analytical
expressions are available. Otherwise a numerical integral is computed.
"""
function integral(f::TimeProfile, t = echotime(f))
    quadgk(f, zero(t), t)
end

function integral(f::PGSE, t = echotime(f))
    δ, Δ = f.δ, f.Δ
    (zero(t) ≤ t < δ) * t + (δ ≤ t < Δ + δ) * δ - (Δ ≤ t < Δ + δ) * (t - Δ)
end

function integral(f::DoublePGSE, t = echotime(f))
    δ, Δ, p = f.δ, f.Δ, f.tpause
    tmid = p + Δ + δ
    if zero(t) ≤ t < δ
        t
    elseif δ ≤ t < Δ
        δ
    elseif Δ ≤ t < Δ + δ
        δ - (t - Δ)
    elseif tmid ≤ t < tmid + δ
        t - tmid
    elseif tmid + δ ≤ t < tmid + Δ
        δ
    elseif tmid + Δ ≤ t < tmid + Δ + δ
        δ - (t - Δ - tmid)
    else
        zero(t)
    end
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


"""
    b = bvalue_no_q(f)

Compute the time profile contribution to the b-value. To obtain the b-value,
multiply the result by `q^2 = (γg)^2`, where `γ` is the gyromagnetic ratio of the
water proton, `g` is the gradient amplitude, and `b = q^2 * bvalue_no_q(f)`.
"""
function bvalue_no_q(f::TimeProfile{T}) where {T}
    quadgk(s -> integral(f, s)^2, zero(T), echotime(f))
end

function bvalue_no_q(f::PGSE)
    f.δ^2 * (f.Δ - f.δ / 3)
end

function bvalue_no_q(f::DoublePGSE)
    2f.δ^2 * (f.Δ - f.δ / 3)
end

function bvalue_no_q(f::CosOGSE)
    f.δ^3 / (2π * f.nperiod)^2
end

function bvalue_no_q(f::SinOGSE)
    3f.δ^3 / (2π * f.nperiod)^2
end


"""
    intervals(f)

Get characteristic intervals of the time profile `f`.
"""
function intervals(f::TimeProfile)
    error("Not implemented")
end

function intervals(f::Union{PGSE,CosOGSE,SinOGSE})
    δ, Δ = f.δ, f.Δ
    if δ < Δ
        [zero(δ), δ, Δ, Δ + δ]
    else
        # No pause between pulses
        [zero(δ), δ, 2f.δ]
    end
end

function intervals(f::DoublePGSE)
    δ, Δ, p = f.δ, f.Δ, f.tpause
    if δ < Δ
        i = [zero(δ), δ, Δ, Δ + δ]
    else
        # No pause between pulses
        i = [zero(δ), δ, 2f.δ, 3f.δ, 4f.δ]
    end
    if p > zero(p)
        [i; i .+ p .+ Δ .+ δ]
    else
        [i; i[2:end] .+ p .+ Δ .+ δ]
    end
end


"""
    TE = echotime(f)

Get echo time `TE` of the time profile `f`, which is the end of the last characteristic
interval.
"""
function echotime(f::TimeProfile)
    intervals(f)[end]
end


"""
    get_interval(f, t)

Get the characteristic interval number of time profile `f` containing time `t`.
If `intervals[i] ≤ t < intervals[i+1]`, `i` is returned.
"""
function get_interval(f::TimeProfile, t)
    i = intervals(f)
    findfirst(i[1:end-1] .≤ t .< i[2:end])
end


"""
    markers = constant_intervals(f)

Get markers for constant intervals. Default is false on all intervals.
"""
function constant_intervals(f::TimeProfile)
    fill(false, length(intervals(f)) - 1)
end

function constant_intervals(f::Union{PGSE,DoublePGSE})
    fill(true, length(intervals(f)) - 1)
end

function constant_intervals(f::Union{CosOGSE,SinOGSE})
    if f.δ < f.Δ
        # OGSE is constant between the pulses
        c = [false, true, false]
    else
        # No pause between pulses
        c = [false, false]
    end
    c
end


"""
    is_constant(f, t)

Return `true` if time profile `f` is constant on the characteristic interval containing or
starting at `t`.
"""
function is_constant(f::TimeProfile, t)
    constant_intervals(f)[get_interval(f, t)]
end
