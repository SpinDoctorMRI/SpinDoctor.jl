# Floating point type for simulations
T = Float64

# Name for saving meshfiles and data
name = "disks"
meshdir = joinpath("meshfiles", name)

# Geometrical setup
setup = DiskSetup{T}(;
    ncell = 3,
    nsidewall = 30,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    # layersizes = [0.6, 1.0],
    ecs = ConvexHullECS{T}(; margin = 2.0),
    # ecs = BoxECS{T}(0.5),
    refinement = 0.1,
)

# Get compartimentalized coefficient vectors
nlayer = length(setup.layersizes)
coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(2) for _ = 1:nlayer], ecs = 0.002 * I(2)),
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
