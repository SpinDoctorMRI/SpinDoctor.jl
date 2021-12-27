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
