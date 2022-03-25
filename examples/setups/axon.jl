# Floating point type for simulations
T = Float64

# Geometrical setup
setup = CylinderSetup{T}(;
    name = "cylinders/axon",
    ncell = 1,
    r_range = (; rmin = 5.0, rmax =  5.0),
    d_range = (; dmin = 0.2, dmax =  0.3),
    deform_angle = (; bend = 0.0, twist =  0.0),
    height = 50.0,
    include_in = true,
    in_ratio = 0.6,
    ecs_shape = :convex_hull,
    ecs_ratio = 0.5,
    refinement = 1,
    D = (; in = 2e-3 * I(3), out = 2e-3 * I(3), ecs = 2e-3 * I(3)),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    T₂ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    κ = SA[1e-4, 1e-4, 0.0, 0.0, 0.0],
)