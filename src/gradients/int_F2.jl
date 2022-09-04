"""
    int_F²(f)

Compute the time profile contribution to the ``b``-value given the time profile `f`. The
``b``-value is given by

```julia
b = γ^2 * g^2 * int_F²(f)
```

where `γ` is the gyromagnetic ratio of the water proton and `g` is the gradient amplitude.
"""
int_F²(f::AbstractTimeProfile{T}) where {T} =
    quadgk(s -> integral(f, s)^2, zero(T), echotime(f))
int_F²(f::PGSE) = f.δ^2 * (f.Δ - f.δ / 3)
int_F²(f::DoublePGSE) = 2f.δ^2 * (f.Δ - f.δ / 3)
int_F²(f::CosOGSE) = f.δ^3 / (2π * f.nperiod)^2
int_F²(f::SinOGSE) = 3f.δ^3 / (2π * f.nperiod)^2
