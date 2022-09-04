"""
    AbstractTimeProfile

General time profile type (gradient sequence).

A `AbstractTimeProfile` should implement a calling method returning the value of the time profile at
time `t`. It can optionally overwrite the [`echotime`](@ref), [`integral`](@ref), and
[`int_F²`](@ref) methods, the first providing a default value and the latter two computing
numerical approximation to the integral quantities if not overwritten.
"""
abstract type AbstractTimeProfile{T} <: Function end

"""
    f = PGSE(δ, Δ)

Pulsed Gradient Spin Echo sequence with pulse duration `δ` and time between
pulses `Δ-δ`.
"""
struct PGSE{T} <: AbstractTimeProfile{T}
    δ::T
    Δ::T
    function PGSE(δ::T, Δ::T) where {T}
        δ ≤ Δ || error("δ > Δ")
        new{T}(δ, Δ)
    end
end

"""
    f = DoublePGSE(δ, Δ)

Double Pulsed Gradient Spin Echo sequence with four pulses of duration `δ`,
separated by pauses of duration `Δ-δ`, `0`, and `Δ-δ` repsectively.
"""
struct DoublePGSE{T} <: AbstractTimeProfile{T}
    δ::T
    Δ::T
    tpause::T
    function DoublePGSE(δ::T, Δ::T, tpause::T) where {T}
        δ ≤ Δ || error("δ > Δ")
        tpause ≥ 0 || error("tpause < 0")
        new{T}(δ, Δ, tpause)
    end
end

"""
    f = CosOGSE(δ, Δ, nperiod)

Oscillating Gradient Spin Echo sequence with two cos-pulses of duration `δ`
separated by a pause of duration `Δ-δ` for `nperiod` periods per pulse.
"""
struct CosOGSE{T} <: AbstractTimeProfile{T}
    δ::T
    Δ::T
    nperiod::Int
    function CosOGSE(δ::T, Δ::T, n) where {T}
        δ ≤ Δ || error("δ > Δ")
        n > 0 || error("nperiod ≤ 0")
        new{T}(δ, Δ, n)
    end
end

"""
    f = SinOGSE(δ, Δ, nperiod)

Oscillating Gradient Spin Echo sequence with two sin-pulses of duration `δ`
separated by a pause of duration `Δ-δ` for `nperiod` periods per pulse.
"""
struct SinOGSE{T} <: AbstractTimeProfile{T}
    δ::T
    Δ::T
    nperiod::Int
    function SinOGSE(δ::T, Δ::T, n) where {T}
        δ ≤ Δ || error("δ > Δ")
        n > 0 || error("nperiod ≤ 0")
        new{T}(δ, Δ, n)
    end
end

"""
    ft = (f::TimeProfile)(t)

Return the time profile value at time `t`.
"""
(f::AbstractTimeProfile)(t) = error("Not implemented")

function (f::PGSE)(t)
    if 0 ≤ t < f.δ
        one(t)
    elseif f.Δ ≤ t < f.Δ + f.δ
        -one(t)
    else
        zero(t)
    end
end

function (f::DoublePGSE)(t)
    δ, Δ, p = f.δ, f.Δ, f.tpause
    if 0 ≤ t < δ
        one(t)
    elseif Δ ≤ t < Δ + δ
        -one(t)
    elseif 0 ≤ t - (p + Δ + δ) < δ
        one(t)
    elseif Δ ≤ t - (p + Δ + δ) < Δ + δ
        -one(t)
    else
        zero(t)
    end
end

function (f::CosOGSE)(t)
    δ, Δ, n = f.δ, f.Δ, f.nperiod
    if 0 ≤ t < δ
        cos(2π * n * t / δ)
    elseif Δ ≤ t < Δ + δ
        -cos(2π * n * (t - Δ) / δ)
    else
        zero(t)
    end
end

function (f::SinOGSE)(t)
    δ, Δ, n = f.δ, f.Δ, f.nperiod
    if 0 ≤ t < δ
        sin(2π * n * t / δ)
    elseif Δ ≤ t < Δ + δ
        -sin(2π * n * (t - Δ) / δ)
    else
        zero(t)
    end
end
