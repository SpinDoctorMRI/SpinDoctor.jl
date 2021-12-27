# Floating point type for simulations
T = Float64

# Geometrical setup
setup = CylinderSetup{T}(;
    name = "cylinders/axon",
    ncell = 1,
    rmin = 5.0,
    rmax = 5.0,
    dmin = 0.2,
    dmax = 0.3,
    bend = 0.0,
    twist = 0.0,
    height = 50.0,
    include_in = true,
    in_ratio = 0.6,
    ecs_shape = "convex_hull",
    ecs_ratio = 0.5,
    refinement = 1,
)

# Get compartimentalized coefficient vectors
coeffs = coefficients(
    setup;
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out= Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)
