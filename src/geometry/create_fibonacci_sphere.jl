function create_fibonacci_sphere(npoint)
    φ = (1 + √5) / 2
    inds = collect(0:npoint-1)'
    θ = 2π / φ * inds
    ϕ = acos.(1 .- 2 .* (inds .+ 0.5) ./ npoint)

    # Create Fibonacci directions (uniformly distributed on unit sphere)
    points = [
        cos.(θ) .* sin.(ϕ)
        sin.(θ) .* sin.(ϕ)
        cos.(ϕ)
    ]
end
