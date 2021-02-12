cellsetup = SpinDoctor.CellSetup(
    name = "slabs",
    shape = "cylinder",
    ncell = 30,
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
    refinement = 0.1
)

domainsetup = SpinDoctor.DomainSetup(
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

experiment = SpinDoctor.ExperimentSetup(
    ndir = 1,
    flat_dirs = false,
    direction = [0.; 0.; 1.],
    sequences = [SpinDoctor.PGSE(2000., 6000.)],
    values = [10000.],
    values_type = 'b',
    btpde = SpinDoctor.BTPDE(odesolver=SpinDoctor.Trapezoid(), nsave=1),
    # mf = MF(length_scale=3, neig_max=400)
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
