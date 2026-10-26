"""
    AbstractTimeProfile

General time profile type (gradient sequence).

A `AbstractTimeProfile` should implement a calling method returning the value of the time
profile at time `t`, and [`intervals`](@ref), returning the characteristic times at which
the profile may be discontinuous (the last one being the echo time). It can optionally
overwrite the [`echotime`](@ref), [`integral`](@ref), and [`int_F²`](@ref) methods, the
first providing a default value and the latter two computing numerical approximation to the
integral quantities if not overwritten.
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

"""
    AbstractGradient{T}

Magnetic field gradient ``\\vec{g}(t) \\in \\mathbb{R}^3``.
"""
abstract type AbstractGradient{T} end

"""
    GeneralGradient(gvec, TE)

General gradient sequence ``\\vec{g}(t) \\in \\mathbb{R}^3``. The direction and amplitude may
vary in time.
"""
Base.@kwdef struct GeneralGradient{T,F}
    gvec::F
    TE::T
end

"""
    ScalarGradient(dir, profile, amplitude)

Gradient sequence with amplitude `amplitude`, normal direction `dir` and scalar time profile
`profile`.

The direction is constant, while the amplitude is controlled by the time profile.

The direction is normalized upon construction.
"""
struct ScalarGradient{T,P<:AbstractTimeProfile{T}}
    dir::Vector{T}
    profile::P
    amplitude::T
    ScalarGradient(dir, profile, amplitude) =
        new{typeof(amplitude),typeof(profile)}(dir / norm(dir), profile, amplitude)
end

(g::GeneralGradient)(t) = g.gvec(t)
(g::ScalarGradient)(t) = g.amplitude * g.profile(t) * g.dir

"""
    intervals(f::AbstractTimeProfile)

Get characteristic intervals of the time profile `f`.

    intervals(g::AbstractGradient)

Get characteristic intervals of the time profile `f`.
"""
function intervals end

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

intervals(g::ScalarGradient) = intervals(g.profile)
intervals(g::GeneralGradient) = [0, g.TE]

"""
    echotime(f::AbstractTimeProfile)

Get echo time ``TE`` of the time profile `f`, which is the end of the last characteristic
interval.

    echotime(grad::AbstractGradient)

Get echo time ``TE`` of `gradient`.
"""
function echotime end

echotime(f::AbstractTimeProfile) = intervals(f)[end]
echotime(g::GeneralGradient) = g.TE
echotime(g::ScalarGradient) = echotime(g.profile)

"""
    isconstant(profile)

Return `true` if the time profile `profile` is interval-wise constant.
"""
function isconstant end

isconstant(::AbstractTimeProfile) = false
isconstant(::PGSE) = true
isconstant(::DoublePGSE) = true

"""
    integral(f, t = echotime(f))

Integral of time profile `f` between `0` and `t`. Unless specified, the echotime
is used as the upper integral limit.

For the `PGSE`, `SinOGSE`, `CosOGSE` and `DoublePGSE` sequences, analytical
expressions are available. Otherwise a numerical integral is computed.
"""
function integral(f::AbstractTimeProfile, t = echotime(f))
    # Split the integral at the characteristic times, where `f` may jump
    breaks = filter(s -> zero(t) < s < t, intervals(f))
    first(quadgk(f, zero(t), breaks..., t))
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
    int_F²(f)

Compute the time profile contribution to the ``b``-value given the time profile `f`. The
``b``-value is given by

```julia
b = γ^2 * g^2 * int_F²(f)
```

where `γ` is the gyromagnetic ratio of the water proton and `g` is the gradient amplitude.
"""
int_F²(f::AbstractTimeProfile) = first(quadgk(s -> integral(f, s)^2, intervals(f)...))
int_F²(f::PGSE) = f.δ^2 * (f.Δ - f.δ / 3)
int_F²(f::DoublePGSE) = 2f.δ^2 * (f.Δ - f.δ / 3)
int_F²(f::CosOGSE) = f.δ^3 / (2π * f.nperiod)^2
int_F²(f::SinOGSE) = 3f.δ^3 / (2π * f.nperiod)^2

"""
    diffusion_time(f)

Diffusion time of time profile `f`.
"""
function diffusion_time end

diffusion_time(f::PGSE) = f.Δ - f.δ / 3
diffusion_time(f::DoublePGSE) = 2(f.Δ - f.δ / 3)
diffusion_time(f::CosOGSE) = 1 / 8 * f.δ / f.nperiod
diffusion_time(f::SinOGSE) = 3 / 8 * f.δ / f.nperiod
