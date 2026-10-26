# Gradient directions

"""
    unitcircle(ndirection; half = false, normal = [0, 0, 1])

Create `ndirection` directions (columns) uniformly distributed on the unit circle in the
plane orthogonal to `normal`.

If `half` is `true`, the directions lie on a half-circle instead of the whole circle.
"""
function unitcircle(ndirection; half = false, normal = [0, 0, 1])
    angles = 2pi * collect(0:(ndirection-1))' / ((1 + half) * ndirection)
    directions = [cos.(angles); sin.(angles); zeros(1, ndirection)]
    reflect_z_to(normal) * directions
end

"""
    unitsphere(ndirection; half = false, normal = [0, 0, 1])

Create `ndirection` directions (columns) uniformly distributed on the unit sphere
(Fibonacci lattice).

If `half` is `true`, the directions are distributed on the hemisphere defined by `normal`
instead.
"""
function unitsphere(ndirection; half = false, normal = [0, 0, 1])
    if half
        reflect_z_to(normal) * fibonacci_sphere(ndirection; half = true)
    else
        fibonacci_sphere(ndirection)
    end
end

"""
    fibonacci_sphere(npoint; half = false)

Create `npoint` evenly distributed points on the unit sphere (or on the upper hemisphere
`z > 0` if `half` is `true`).
"""
function fibonacci_sphere(npoint; half = false)
    golden = (1 + sqrt(5)) / 2
    inds = collect(0:(npoint-1))'
    azimuth = 2pi / golden * inds
    polar = acos.(1 .- (2 - half) .* (inds .+ 0.5) ./ npoint)
    [
        cos.(azimuth) .* sin.(polar)
        sin.(azimuth) .* sin.(polar)
        cos.(polar)
    ]
end

# Reflection mapping e_z to the unit vector along `normal` (the identity if `normal` is e_z)
function reflect_z_to(normal)
    n = normal / norm(normal)
    isapprox(n, [0, 0, 1]) && return Matrix(1.0I, 3, 3)
    c = n + [0, 0, 1]
    norm(c) > sqrt(eps()) || return -Matrix(1.0I, 3, 3)
    2 * (c * c') / (c' * c) - I
end
