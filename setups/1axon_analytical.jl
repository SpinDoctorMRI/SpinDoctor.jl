setup = CylinderSetup(
    name = "cylinders/1axons_analytical",
    ncell = 1,
    rmin = 5.0,
    rmax = 5.0,
    dmin = 0.2,
    dmax = 0.3,
    bend = 0.0,
    twist = 0.0,
    height = 50.0,
    include_in = true,
    in_ratio = 0.6,
    ecs_shape = "convex_hull",
    ecs_ratio = 0.5,
    refinement = 1,
    D_in = 0.002 * I(3),
    D_out = 0.002 * I(3),
    D_ecs = 0.002 * I(3),
    T₂_in = Inf,
    T₂_out = Inf,
    T₂_ecs = Inf,
    ρ_in = 1.0,
    ρ_out = 1.0,
    ρ_ecs = 1.0,
    κ_in_out = 1e-4,
    κ_out_ecs = 1e-4,
    κ_in = 0.0,
    κ_out = 0.0,
    κ_ecs = 0.0,
)

experiment = Experiment(
    gradient = (
        directions = create_directions([1.0; 0.0; 0.0]),
        sequences = [PGSE(5000, 10000), PGSE(10000, 100000)],
        values = 0.0:500:10000,
        values_type = "b",
    ),
    btpde = (
        odesolver = Rodas5(autodiff = :false),
        reltol = 1e-4,
        abstol = 1e-6,
        nsave = 1,
    ),
    mf = (length_scale = 3, neig_max = 400, ninterval = 500),
    analytical = (length_scale = 0.3, eigstep = 1e-8),
)


# ImplicitEuler()
# Trapezoid()
# ABDF2()
# QNDF()
# QNDF1()
# QNDF2()
# QBDF()
# QBDF1()
# QBDF2()
# ROS3P()
# Rosenbrock23(autodiff = :false)
# ROS34PW1a()
# ROS34PW3()
# Rodas4(autodiff = :false)
# Rodas4P(autodiff = :false)
# Rodas5(autodiff = :false)
# Kvaerno3()
# KenCarp4()
# Cash4()
# Tsit5()
# Vern7()
# VCABM()
# BS3()
# DP5()
# DP8()