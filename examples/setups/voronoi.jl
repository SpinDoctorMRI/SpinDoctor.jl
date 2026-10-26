# Name for saving meshfiles and data
name = "voronoi"
savedir = joinpath("data", name)

# Geometrical setup: polyhedral cells filling 85 % of a periodic box, cut by the box faces.
# The mesh is periodic: use `assemble_matrices(model; periodic = true)` for periodic
# boundary conditions.
setup = VoronoiSetup(; ncell = 10, boxsize = 15.0, fraction = 0.85)

# Get compartimentalized coefficient vectors
coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3)], ecs = 0.002 * I(3)),
    T2 = (; cell = [Inf], ecs = Inf),
    rho = (; cell = [1.0], ecs = 1.0),
    kappa = (;
        cell_interfaces = Float64[],
        cell_boundaries = [0],
        cell_ecs = 1e-4,
        ecs = 0,
    ),
    gamma = 2.67513e-4,
)
