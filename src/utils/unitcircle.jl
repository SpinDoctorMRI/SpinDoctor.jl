"""
    unitcircle(ndirection; half = false, normal = [0, 0, 1])

Create `ndirection` directions unformly distributed on the unit circle defined by `normal`.

If `half` is `true`, the directions lie on a half-circle instead of the whole circle. 
"""
function unitcircle(n; half = false, normal = [0, 0, 1])
    # Create directions (unformly distributed on unit circle)
    θ = 2π * collect(0:n-1)' / (1 + half)n
    directions = [
        cos.(θ)
        sin.(θ)
        zeros(1, n)
    ]

    # Create rotation matrix to transform normal to e_z
    normal = normal / norm(normal)
    c = normal + [0, 0, 1]
    R = 2(c * c') / c'c - I(3)
    normal ≈ [0, 0, 1] || (directions = R'directions)

    directions
end
