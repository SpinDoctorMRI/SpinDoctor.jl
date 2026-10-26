# Name for saving meshfiles and data
name = "packed_cylinders"
savedir = joinpath("data", name)

# Geometrical setup: a dense packing of myelinated axons in a periodic box (70 % of the
# area), cut by the box faces. The mesh is periodic: use
# `assemble_matrices(model; periodic = true)` for periodic boundary conditions.
setup = CylinderSetup(;
    ncell = 20,
    nsidewall = 16,
    rmin = 1.0,
    rmax = 3.0,
    dmin = 0.1,
    layersizes = [0.7, 1.0],
    ecs = PeriodicBoxECS(; fraction = 0.7),
    height = 5.0,
)

# Get compartimentalized coefficient vectors
nlayer = length(setup.groundsetup.layersizes)
coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3), 0.0005 * I(3)], ecs = 0.002 * I(3)),
    T2 = (; cell = fill(Inf, nlayer), ecs = Inf),
    rho = (; cell = fill(1.0, nlayer), ecs = 1.0),
    kappa = (;
        cell_interfaces = fill(1e-5, nlayer - 1),
        cell_boundaries = fill(0, nlayer),
        cell_ecs = 1e-5,
        ecs = 0,
    ),
    gamma = 2.67513e-4,
)
