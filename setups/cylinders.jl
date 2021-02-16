cellsetup = CellSetup(
    name = "slabs",
    shape = "cylinder",
    ncell = 5,
    rmin = 1.,
    rmax = 8.,
    dmin = 0.2,
    dmax = 0.3,
    deformation = (0., 0.), #(0.05, 2π/3), #
    height = 1.,
    include_nucleus = true,
    nucleus_radiusratio = 0.7,
    include_ecs = true,
    ecs_shape = "convexhull",
    ecs_gap = 0.3,
    refinement = 0.1,
)

domainsetup = DomainSetup(
    diffusivity_in = 0.002,
    diffusivity_out = 0.002,
    diffusivity_ecs = 0.002,
    relaxation_in = Inf,
    relaxation_out = Inf,
    relaxation_ecs = Inf,
    initial_density_in = 1.,
    initial_density_out = 1.,
    initial_density_ecs = 1.,
    permeability_in_out = 1e-4,
    permeability_out_ecs = 1e-4,
)

experiment = ExperimentSetup(
    ndirection = 1,
    flat_dirs = false,
    direction = [1.; 0.; 0.],
    sequences = [PGSE(2000., 6000.)],
    values = [50., 100., 200., 1000., 4000., 10000.],
    values_type = 'b',
    btpde = BTPDE(odesolver=Trapezoid(), nsave=1),
    # mf = MF(length_scale=3, neig_max=400),
)

# ImplicitEuler
# Trapezoid
# ABDF2
# QNDF
# QNDF1
# QNDF2
# QBDF
# QBDF1
# QBDF2
