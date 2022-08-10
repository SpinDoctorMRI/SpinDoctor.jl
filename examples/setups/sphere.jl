# Floating point type for simulations
T = Float64

# Name for saving meshfiles and data
name = "sphere"
meshdir = joinpath("meshfiles", name)

# Geometrical setup
setup = SphereSetup{T}(;
    ncell = 1,
    layersizes = [0.6, 1.0],
    rmin = 5.0,
    rmax = 5.0,
    dmin = 0.2,
    dmax = 0.3,
    nsidewall = 200,
    ecs = NoECS(),
    refinement = 0.5,
)

# Get compartimentalized coefficient vectors
nlayer = length(setup.layersizes)
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
