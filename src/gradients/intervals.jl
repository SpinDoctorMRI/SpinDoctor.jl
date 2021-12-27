"""
    intervals(f)

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
