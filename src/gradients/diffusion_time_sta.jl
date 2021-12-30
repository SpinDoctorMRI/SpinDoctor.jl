"""
    diffusion_time_sta(f)

Get the short term approximation of the diffusion time of the time profile `f`.
"""
function diffusion_time_sta end

function diffusion_time_sta(f::PGSE)
    (; δ, Δ) = f
    out =
        4 / 35 * ((Δ + δ)^(7 / 2) + (Δ - δ)^(7 / 2) - 2Δ^(7 / 2) - 2δ^(7 / 2)) /
        (δ^2 * (Δ - δ / 3))
    out^2
end

function diffusion_time_sta(f::DoublePGSE)
    (; δ, Δ, tpause) = f
    tm = δ + tpause
    out =
        2 / 35 * (
            (2Δ + tm + δ)^(7 / 2) +
            (2Δ + tm - δ)^(7 / 2) +
            (tm + δ)^(7 / 2) +
            (tm - δ)^(7 / 2) - 2(2Δ + tm)^(7 / 2) - 2(Δ + tm + δ)^(7 / 2) -
            2(Δ + tm - δ)^(7 / 2) +
            2(Δ + δ)^(7 / 2) +
            2(Δ - δ)^(7 / 2) - 2tm^(7 / 2) + 4(Δ + tm)^(7 / 2) - 4Δ^(7 / 2) - 4δ^(7 / 2)
        ) / (δ^2 * (Δ - δ / 3))
    out^2
end

function diffusion_time_sta(f::CosOGSE)
    (; δ, Δ) = f
    n = f.nperiod
    S = sin(2π * n * Δ / δ)
    C = cos(2π * n * Δ / δ)
    FS1 = fresnels(2√n)
    FS2 = fresnels(2√(Δ * n / δ))
    FS3 = fresnels(2√((Δ + δ) * n / δ))
    FS4 = fresnels(2√((Δ - δ) * n / δ))
    FC1 = fresnelc(2√n)
    FC2 = fresnelc(2√(Δ * n / δ))
    FC3 = fresnelc(2√((Δ + δ) * n / δ))
    FC4 = fresnelc(2√((Δ - δ) * n / δ))
    out =
        3 / 8 / √(n * δ) * (
            2Δ * FC2 * C + δ * FC4 * C - δ * FC3 * C - Δ * FC3 * C - Δ * FC4 * C +
            2Δ * FS2 * S +
            δ * FS4 * S - δ * FS3 * S - Δ * FS3 * S - Δ * FS4 * S + 2δ * FC1
        ) +
        9√δ / 32π / n^(3 / 2) *
        (+2FS2 * C - FS3 * C - FS4 * C - 2FC2 * S + FC3 * S + FC4 * S + 2FS1)
    out^2
end

function diffusion_time_sta(f::SinOGSE)
    (; δ, Δ) = f
    n = f.nperiod
    S = sin(2π * n * Δ / δ)
    C = cos(2π * n * Δ / δ)
    FS1 = fresnels(2√n)
    FS2 = fresnels(2√(Δ * n / δ))
    FS3 = fresnels(2√((Δ + δ) * n / δ))
    FS4 = fresnels(2√((Δ - δ) * n / δ))
    FC1 = fresnelc(2√n)
    FC2 = fresnelc(2√(Δ * n / δ))
    FC3 = fresnelc(2√((Δ + δ) * n / δ))
    FC4 = fresnelc(2√((Δ - δ) * n / δ))
    out =
        (
            64π * n^(3 / 2) * Δ^(3 / 2) +
            64π * n^(3 / 2) * δ^(3 / 2) +
            42δ^(3 / 2) * FS2 * C +
            42δ^(3 / 2) * FS1 - 42δ^(3 / 2) * FC2 * S + 32π * n^(3 / 2) * √(Δ - δ) * δ -
            32π * n^(3 / 2) * √(Δ + δ) * Δ - 32π * n^(3 / 2) * √(Δ + δ) * δ -
            32π * n^(3 / 2) * √(Δ - δ) * Δ +
            24π * n * Δ * √δ * FC2 * C +
            24π * n * Δ * √δ * FS2 * S +
            24π * n * δ^(3 / 2) * FC1 +
            21δ^(3 / 2) * FC3 * S +
            21δ^(3 / 2) * FC4 * S - 21δ^(3 / 2) * FS3 * C - 21δ^(3 / 2) * FS4 * C +
            12π * n * δ^(3 / 2) * FS4 * S +
            12π * n * δ^(3 / 2) * FC4 * C - 12π * n * δ^(3 / 2) * FC3 * C -
            12π * n * δ^(3 / 2) * FS3 * S - 12π * n * Δ * √δ * FC4 * C -
            12π * n * Δ * √δ * FC3 * C - 12π * n * Δ * √δ * FS3 * S -
            12π * n * Δ * √δ * FS4 * S
        ) / (96π * n^(3 / 2) * δ)
    out^2
end
