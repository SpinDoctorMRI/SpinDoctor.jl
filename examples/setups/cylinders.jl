# Floating point type for simulations
T = Float64

# Geometrical setup
setup = CylinderSetup{T}(;
    name = "cylinders/somecylinders",
    ncell = 10,
    r_range = (; rmin = 2.0, rmax =  6.0),
    d_range = (; dmin = 0.2, dmax =  0.3),
    height = 20.0,
    deform_angle = (; bend = 0.0, twist =  0.0),
    include_in = false,
    in_ratio = 0.6,
    ecs_shape = :convex_hull,
    ecs_ratio = 0.5,
    refinement = 0.5,
    D = (; in = 2e-3 * I(3), out = 2e-3 * I(3), ecs = 2e-3 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 0e-4, out_ecs = 0e-4, in = 0.0, out = 0.0, ecs = 0.0),
)