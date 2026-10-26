# Time profiles and magnetic field gradient sequences

"""
    AbstractTimeProfile{T}

Time profile `f(t)` of a gradient sequence (normalized to a maximum of one).

A time profile implements a calling method returning its value at time `t`, and
[`intervals`](@ref), returning the characteristic times at which the profile may be
discontinuous (the last one being the echo time). It can optionally implement
[`integral`](@ref) and [`int_F2`](@ref), which are otherwise computed numerically.
"""
abstract type AbstractTimeProfile{T} <: Function end

"""
    PGSE(delta, Delta)

Pulsed gradient spin echo (PGSE) time profile with pulse duration `delta` and time
`Delta` between the starts of the two pulses.
"""
struct PGSE{T} <: AbstractTimeProfile{T}
    delta::T
    Delta::T
    function PGSE(delta::T, Delta::T) where {T<:AbstractFloat}
        delta <= Delta || error("delta > Delta")
        new{T}(delta, Delta)
    end
end
PGSE(delta, Delta) = PGSE(promote(float(delta), float(Delta))...)

"""
    DoublePGSE(delta, Delta, tpause)

Double pulsed gradient spin echo time profile with four pulses of duration `delta`,
separated by pauses of duration `Delta - delta`, `tpause`, and `Delta - delta`
respectively.
"""
struct DoublePGSE{T} <: AbstractTimeProfile{T}
    delta::T
    Delta::T
    tpause::T
    function DoublePGSE(delta::T, Delta::T, tpause::T) where {T<:AbstractFloat}
        delta <= Delta || error("delta > Delta")
        tpause >= 0 || error("tpause < 0")
        new{T}(delta, Delta, tpause)
    end
end
DoublePGSE(delta, Delta, tpause) =
    DoublePGSE(promote(float(delta), float(Delta), float(tpause))...)

"""
    CosOGSE(delta, Delta, nperiod)

Oscillating gradient spin echo time profile with two cosine pulses of duration `delta`,
starting at times `0` and `Delta`, with `nperiod` periods per pulse.
"""
struct CosOGSE{T} <: AbstractTimeProfile{T}
    delta::T
    Delta::T
    nperiod::Int
    function CosOGSE(delta::T, Delta::T, nperiod) where {T<:AbstractFloat}
        delta <= Delta || error("delta > Delta")
        nperiod > 0 || error("nperiod <= 0")
        new{T}(delta, Delta, nperiod)
    end
end
CosOGSE(delta, Delta, nperiod) = CosOGSE(promote(float(delta), float(Delta))..., nperiod)

"""
    SinOGSE(delta, Delta, nperiod)

Oscillating gradient spin echo time profile with two sine pulses of duration `delta`,
starting at times `0` and `Delta`, with `nperiod` periods per pulse.
"""
struct SinOGSE{T} <: AbstractTimeProfile{T}
    delta::T
    Delta::T
    nperiod::Int
    function SinOGSE(delta::T, Delta::T, nperiod) where {T<:AbstractFloat}
        delta <= Delta || error("delta > Delta")
        nperiod > 0 || error("nperiod <= 0")
        new{T}(delta, Delta, nperiod)
    end
end
SinOGSE(delta, Delta, nperiod) = SinOGSE(promote(float(delta), float(Delta))..., nperiod)

(f::AbstractTimeProfile)(t) = error("Not implemented")

function (f::PGSE)(t)
    (; delta, Delta) = f
    if 0 <= t < delta
        one(t)
    elseif Delta <= t < Delta + delta
        -one(t)
    else
        zero(t)
    end
end

function (f::DoublePGSE)(t)
    (; delta, Delta, tpause) = f
    tmid = tpause + Delta + delta
    if 0 <= t < delta
        one(t)
    elseif Delta <= t < Delta + delta
        -one(t)
    elseif 0 <= t - tmid < delta
        one(t)
    elseif Delta <= t - tmid < Delta + delta
        -one(t)
    else
        zero(t)
    end
end

function (f::CosOGSE)(t)
    (; delta, Delta, nperiod) = f
    if 0 <= t < delta
        cos(2pi * nperiod * t / delta)
    elseif Delta <= t < Delta + delta
        -cos(2pi * nperiod * (t - Delta) / delta)
    else
        zero(t)
    end
end

function (f::SinOGSE)(t)
    (; delta, Delta, nperiod) = f
    if 0 <= t < delta
        sin(2pi * nperiod * t / delta)
    elseif Delta <= t < Delta + delta
        -sin(2pi * nperiod * (t - Delta) / delta)
    else
        zero(t)
    end
end

"""
    intervals(f)

Characteristic times of the time profile (or gradient) `f`: the profile is smooth between
two consecutive times. The first time is zero and the last one is the echo time.
"""
function intervals end

function intervals(f::Union{PGSE,CosOGSE,SinOGSE})
    (; delta, Delta) = f
    if delta < Delta
        [zero(delta), delta, Delta, Delta + delta]
    else
        # No pause between pulses
        [zero(delta), delta, 2delta]
    end
end

function intervals(f::DoublePGSE)
    (; delta, Delta, tpause) = f
    if delta < Delta
        i = [zero(delta), delta, Delta, Delta + delta]
    else
        # No pause between pulses
        i = [zero(delta), delta, 2delta, 3delta, 4delta]
    end
    if tpause > 0
        [i; i .+ tpause .+ Delta .+ delta]
    else
        [i; i[2:end] .+ Delta .+ delta]
    end
end

"""
    echotime(f)

Echo time of the time profile (or gradient) `f`, which is the end of the last
characteristic interval.
"""
echotime(f::AbstractTimeProfile) = intervals(f)[end]

"""
    isconstant(f)

Return `true` if the time profile `f` is interval-wise constant.
"""
isconstant(::AbstractTimeProfile) = false
isconstant(::PGSE) = true
isconstant(::DoublePGSE) = true

"""
    integral(f, t = echotime(f))

Integral `F(t)` of the time profile `f` between `0` and `t`.

Closed forms are used for [`PGSE`](@ref), [`DoublePGSE`](@ref), [`CosOGSE`](@ref) and
[`SinOGSE`](@ref). Otherwise the integral is computed numerically.
"""
function integral(f::AbstractTimeProfile, t = echotime(f))
    # Split the integral at the characteristic times, where `f` may jump
    breaks = filter(s -> zero(t) < s < t, intervals(f))
    first(quadgk(f, zero(t), breaks..., t))
end

function integral(f::PGSE, t = echotime(f))
    (; delta, Delta) = f
    (0 <= t < delta) * t + (delta <= t < Delta + delta) * delta -
    (Delta <= t < Delta + delta) * (t - Delta)
end

function integral(f::DoublePGSE, t = echotime(f))
    (; delta, Delta, tpause) = f
    tmid = tpause + Delta + delta
    if 0 <= t < delta
        t
    elseif delta <= t < Delta
        delta
    elseif Delta <= t < Delta + delta
        delta - (t - Delta)
    elseif tmid <= t < tmid + delta
        t - tmid
    elseif tmid + delta <= t < tmid + Delta
        delta
    elseif tmid + Delta <= t < tmid + Delta + delta
        delta - (t - Delta - tmid)
    else
        zero(t)
    end
end

function integral(f::CosOGSE, t = echotime(f))
    (; delta, Delta, nperiod) = f
    a = 2pi * nperiod / delta
    ((t < delta) * sin(a * t) - (Delta <= t) * sin(a * (t - Delta))) / a
end

function integral(f::SinOGSE, t = echotime(f))
    (; delta, Delta, nperiod) = f
    a = 2pi * nperiod / delta
    ((t < delta) * (1 - cos(a * t)) - (Delta <= t) * (1 - cos(a * (t - Delta)))) / a
end

"""
    int_F2(f)

Integral ``\\int_0^{T_E} F(t)^2 \\, \\mathrm{d}t`` of the squared integral `F` of the time
profile `f`.
It relates the gradient amplitude `g` to the b-value: `b = gamma^2 * g^2 * int_F2(f)`.
"""
int_F2(f::AbstractTimeProfile) = first(quadgk(s -> integral(f, s)^2, intervals(f)...))
int_F2(f::PGSE) = f.delta^2 * (f.Delta - f.delta / 3)
int_F2(f::DoublePGSE) = 2 * f.delta^2 * (f.Delta - f.delta / 3)
int_F2(f::CosOGSE) = f.delta^3 / (2pi * f.nperiod)^2
int_F2(f::SinOGSE) = 3 * f.delta^3 / (2pi * f.nperiod)^2

"""
    diffusion_time(f)

Diffusion time of the time profile `f`.
"""
function diffusion_time end

diffusion_time(f::PGSE) = f.Delta - f.delta / 3
diffusion_time(f::DoublePGSE) = 2 * (f.Delta - f.delta / 3)
diffusion_time(f::CosOGSE) = f.delta / f.nperiod / 8
diffusion_time(f::SinOGSE) = 3 * f.delta / f.nperiod / 8

"""
    AbstractGradient{T}

Magnetic field gradient sequence `g(t)` (a vector).
"""
abstract type AbstractGradient{T} end

"""
    ScalarGradient(dir, profile, amplitude)
    ScalarGradient(dir, profile; bvalue, gamma)

Gradient sequence `g(t) = amplitude * profile(t) * dir`, with a constant direction `dir`
(normalized upon construction), a time profile `profile` and an amplitude `amplitude`.

Instead of the amplitude, the b-value may be given together with the gyromagnetic ratio
`gamma`.
"""
struct ScalarGradient{T,P<:AbstractTimeProfile{T}} <: AbstractGradient{T}
    dir::Vector{T}
    profile::P
    amplitude::T
    function ScalarGradient(dir, profile::AbstractTimeProfile{T}, amplitude) where {T}
        new{T,typeof(profile)}(dir / norm(dir), profile, amplitude)
    end
end

ScalarGradient(dir, profile; bvalue, gamma) =
    ScalarGradient(dir, profile, sqrt(bvalue / int_F2(profile)) / gamma)

"""
    GeneralGradient(; gvec, TE)

General gradient sequence `g(t) = gvec(t)`, whose direction and amplitude may vary in time,
with echo time `TE`.
"""
Base.@kwdef struct GeneralGradient{T,F} <: AbstractGradient{T}
    gvec::F
    TE::T
end

(g::ScalarGradient)(t) = g.amplitude * g.profile(t) * g.dir
(g::GeneralGradient)(t) = g.gvec(t)

intervals(g::ScalarGradient) = intervals(g.profile)
intervals(g::GeneralGradient) = [zero(g.TE), g.TE]

echotime(g::ScalarGradient) = echotime(g.profile)
echotime(g::GeneralGradient) = g.TE

"""
    bvalue(gradient, gamma)

b-value of the scalar gradient sequence `gradient` for the gyromagnetic ratio `gamma`.
"""
bvalue(g::ScalarGradient, gamma) = (gamma * g.amplitude)^2 * int_F2(g.profile)
