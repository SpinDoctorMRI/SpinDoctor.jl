# Floating point type for simulations
T = Float64

# Geometrical setup
setup = SphereSetup{T}(
    name = "spheres/somespheres",
    ncell = 4,
    r_range = (; rmin = 3.0, rmax =  5.0),
    d_range = (; dmin = 0.2, dmax =  0.3),
    include_in = false,
    in_ratio = 0.7,
    ecs_shape = :convex_hull,
    ecs_ratio = 0.3,
    refinement = 0.5,
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
)