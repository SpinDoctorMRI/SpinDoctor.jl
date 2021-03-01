cellsetup = CellSetup(
    name = "spindle/whole_neurons/03b_spindle4aACC",
    shape = "cylinder",
    ncell = 1,
    # rmin = 1.,
    # rmax = 8.,
    # dmin = 0.2,
    # dmax = 0.3,
    # deformation = (0., 0.), #(0.05, π/3), #
    # height = 1.,
    include_in = false,
    in_ratio = 0.7,
    include_ecs = false,
    ecs_shape = "convexhull",
    ecs_ratio = 0.3,
    refinement = 0.3
)

domainsetup = DomainSetup(
    σ_in = 0.002,
    σ_out = 0.002,
    σ_ecs = 0.002,
    T₂_in = Inf,
    T₂_out = Inf,
    T₂_ecs = Inf,
    ρ_in = 1.,
    ρ_out = 1.,
    ρ_ecs = 1.,
    κ_in_out = 1e-4,
    κ_out_ecs = 1e-4,
    κ_in = 0.,
    κ_out = 0.,
    κ_ecs = 0.,
)

experiment = ExperimentSetup(
    ndirection = 1,
    flat_dirs = false,
    direction = [1.; 0.; 0.],
    sequences = [PGSE(2000., 6000.)],
    values = [10000.],
    values_type = 'b',
    btpde = BTPDE(odesolver=Trapezoid(), nsave=100),
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
