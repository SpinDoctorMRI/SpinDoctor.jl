# Create setup
setup = Setup(name = "spindle/whole_neurons/03b_spindle4aACC")

# Geometry paramters
setup.geometry = Dict(
    :cell_shape => "neuron",
    :ncell => 1,
    :include_in => false,
    :in_ratio => 0.7,
    :ecs_shape => "no_ecs",
    :ecs_ratio => 0.3,
    :deformation => (0., 0.),
    # :refinement => 0.5,
)

# PDE paramaters
setup.pde = Dict(
    :σ_in => 0.002,
    :σ_out => 0.002,
    :σ_ecs => 0.002,
    :T₂_in => Inf,
    :T₂_out => Inf,
    :T₂_ecs => Inf,
    :ρ_in => 1.,
    :ρ_out => 1.,
    :ρ_ecs => 1.,
    :κ_in_out => 1e-4,
    :κ_out_ecs => 1e-4,
    :κ_in => 0.,
    :κ_out => 0.,
    :κ_ecs => 0.,
)

# Gradient sequences
setup.gradient = Dict(
    :directions => create_directions([1.; 0.; 0.]),
    #:directions => create_directions(10, flat=false),
    :sequences => [PGSE(2000., 6000.)],
    # :values => 0:500:2000,
    :values => [1000.],
    :values_type => 'b',
)

# BTPDE solver
setup.btpde = Dict(
    :odesolver => QBDF(),
    :reltol => 1e-4,
    :abstol => 1e-6,
    :nsave => 1,
)

# MF solver
setup.mf = Dict(
    :length_scale => 3,
    :neig_max => 400,
    :ninterval => 500,
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