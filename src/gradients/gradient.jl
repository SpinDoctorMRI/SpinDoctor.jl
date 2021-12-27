"""
    AbstractGradient{T}

Magnetic field gradient ``\\vec{g}(t)``.
"""
abstract type AbstractGradient{T} end

"""
    GeneralGradient(g⃗)

General gradient sequence `g⃗(t) ∈ R³`. The direcion and amplitude may vary in time.
"""
Base.@kwdef struct GeneralGradient{T,F}
    g⃗::F
    TE::T
end

"""
    ScalarGradient(dir, profile, amplitude)

Gradient sequence with amplitude `amplitude`, normal direction `dir` and scalar time profile
`profile`.

The direction is constant, while the amplitude is controlled by the time profile.

The direction is normalized upon construction.
"""
struct ScalarGradient{T}
    dir::Vector{T}
    profile::TimeProfile{T}
    amplitude::T 
    ScalarGradient(dir, profile, amplitude) =
        new{typeof(amplitude)}(dir / norm(dir), profile, amplitude)
end

(g::GeneralGradient)(t) = g.g⃗(t)
(g::ScalarGradient)(t) = g.amplitude * g.profile(t) * g.dir
