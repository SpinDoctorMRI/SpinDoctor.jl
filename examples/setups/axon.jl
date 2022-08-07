# Floating point type for simulations
T = Float64

# Name for saving meshfiles and data
name = "axon"
meshdir = joinpath("meshfiles", name)

# Geometrical setup
setup = CylinderSetup{T}(;
    ncell = 1,
    nsidewall = 30,
    rmin = 5.0,
    rmax = 5.0,
    layersizes = [0.25, 0.5, 0.75, 1.0],
    ecs_shape = :no_ecs,
    height = 50.0,
    bend = 0.0,
    twist = 0.0,
    refinement = 1.0,
)

# Get compartimentalized coefficient vectors
nlayer = length(groundsetup.layersizes)
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
