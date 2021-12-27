"""
    unitsphere(ndirection; half = false, normal = [0, 0, 1])

Create `ndirection` directions unformly distributed on the unit sphere.

If `half` is `true`, the points will be discributed on the hemisphere defined by `normal`.
"""
function unitsphere(ndirection; half = false, normal = [0, 0, 1])
    # Create Fibonacci distributed directions
    directions = create_fibonacci_sphere(ndirection)
end
