# Floating point type for simulations
T = Float64

# Geometrical setup
groundsetup = DiskSetup{T}(;
    ncell = 3,
    nsidewall = 12,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    layersizes = [1.0],
    ecs_shape = :convex_hull,
    ecs_ratio = 0.5,
    # refinement = 1.0,
)
setup = ExtrusionSetup(;
    name = "cylinders/somecylinders",
    groundsetup,
    height = 5.0,
    bend = 0.0,
    twist = 0.0,
    # refinement = 0.5,
)

# Get compartimentalized coefficient vectors
nlayer = length(setup.groundsetup.layersizes)
coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3) for _ = 1:nlayer], ecs = 0.002 * I(3)),
    T₂ = (; cell = fill(Inf, nlayer), ecs = Inf),
    ρ = (; cell = fill(1.0, nlayer), ecs = 1.0),
    κ = (;
        cell_interfaces = fill(1e-4, nlayer - 1),
        cell_boundaries = fill(0, nlayer),
        cell_ecs = 1e-4,
        ecs = 0,
    ),
    γ = 2.67513e-4,
)
