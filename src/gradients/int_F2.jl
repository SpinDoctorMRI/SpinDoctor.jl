"""
    b = int_F²(f)

Compute the time profile contribution to the b-value. To obtain the b-value,
multiply the result by `q^2 = (γg)^2`, where `γ` is the gyromagnetic ratio of the
water proton, `g` is the gradient amplitude, and `b = q^2 * int_F²(f)`.
"""
function int_F²(f::TimeProfile{T}) where {T}
    quadgk(s -> integral(f, s)^2, zero(T), echotime(f))
end

function int_F²(f::PGSE)
    f.δ^2 * (f.Δ - f.δ / 3)
end

function int_F²(f::DoublePGSE)
    2f.δ^2 * (f.Δ - f.δ / 3)
end

function int_F²(f::CosOGSE)
    f.δ^3 / (2π * f.nperiod)^2
end

function int_F²(f::SinOGSE)
    3f.δ^3 / (2π * f.nperiod)^2
end

