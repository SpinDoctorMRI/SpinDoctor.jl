"""
    create_directions(ndirection, flat)

Create `ndirection` directions unformly distributed in sphere.
If `flat` is true, the directions lie in the plane instead of on the sphere.
"""
function create_directions(ndirection::Int; flat::Bool = false)
    if flat
        # Create directions (unformly distributed on unit circle)
        θ = 2π * collect(0:ndirection-1)' / ndirection
        directions = [
            cos.(θ)
            sin.(θ)
            zeros(1, ndirection)
        ]
    else
        # Create Fibonacci distributed directions
        directions = create_fibonacci_sphere(ndirection)
    end
    directions
end
