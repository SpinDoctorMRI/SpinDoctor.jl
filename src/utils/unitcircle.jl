"""
    unitcircle(ndirection; half = false, normal = [0, 0, 1])

Create `ndirection` directions unformly distributed on the unit circle defined by `normal`.

If `half` is `true`, the directions lie on a half-circle instead of the whole circle. 
"""
function unitcircle(n; half = false, normal = [0, 0, 1])
    # Create directions (unformly distributed on unit circle)
    θ = 2π * collect(0:n-1)' / n
    directions = [
        cos.(θ)
        sin.(θ)
        zeros(1, ndirection)
    ]
end
