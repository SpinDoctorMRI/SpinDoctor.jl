# Floating point type for simulations
T = Float64

# Geometrical setup
setup = PlateSetup{T}(;
    name = "plates/someplates",
    width = 1.0,
    depth = 1.0,
    heights = 0:0.1:1.0,
    bend = 0.0,
    twist = 0.0,
    # refinement = 0.5,
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
