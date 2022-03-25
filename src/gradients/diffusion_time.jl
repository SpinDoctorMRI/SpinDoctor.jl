"""
    diffusion_time(f)

Diffusion time of time profile `f`.
"""
function diffusion_time end

diffusion_time(f::PGSE) = f.Δ - f.δ / 3
diffusion_time(f::DoublePGSE) = 2(f.Δ - f.δ / 3)
diffusion_time(f::CosOGSE) = 1 / 8 * f.δ / f.nperiod
diffusion_time(f::SinOGSE) = 3 / 8 * f.δ / f.nperiod
