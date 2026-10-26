"""
    unitcircle(ndirection; half = false, normal = [0, 0, 1])

Create `ndirection` directions uniformly distributed on the unit circle
defined by `normal`.

If `half` is `true`, the directions lie on a half-circle instead of the whole
circle. 
"""
function unitcircle(n; half = false, normal = [0, 0, 1])
    # Create directions (unformly distributed on unit circle)
    θ = 2π * collect(0:(n-1))' / (1 + half)n
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

"""
    unitsphere(ndirection; half = false, normal = [0, 0, 1])

Create `ndirection` directions uniformly distributed on the unit sphere (Fibonacci
lattice).

If `half` is `true`, the directions are distributed on the hemisphere defined by `normal`
instead.
"""
function unitsphere(ndirection; half = false, normal = [0, 0, 1])
    if half
        # Fibonacci lattice on the upper hemisphere (z > 0)
        φ = (1 + √5) / 2
        inds = collect(0:(ndirection-1))'
        θ = 2π / φ * inds
        ϕ = acos.(1 .- (inds .+ 0.5) ./ ndirection)
        directions = [
            cos.(θ) .* sin.(ϕ)
            sin.(θ) .* sin.(ϕ)
            cos.(ϕ)
        ]

        # Reflect e_z onto `normal`
        normal = normal / norm(normal)
        c = normal + [0, 0, 1]
        if norm(c) > √eps()
            R = 2(c * c') / c'c - I(3)
            directions = R * directions
        else
            directions = -directions
        end
    else
        directions = create_fibonacci_sphere(ndirection)
    end
    directions
end

"""
    create_fibonacci_sphere(npoint)

Create `npoint` evenly distributed points on the sphere.
"""
function create_fibonacci_sphere(npoint)
    φ = (1 + √5) / 2
    inds = collect(0:(npoint-1))'
    θ = 2π / φ * inds
    ϕ = acos.(1 .- 2 .* (inds .+ 0.5) ./ npoint)

    # Create Fibonacci directions (uniformly distributed on unit sphere)
    [
        cos.(θ) .* sin.(ϕ)
        sin.(θ) .* sin.(ϕ)
        cos.(ϕ)
    ]
end
