# Create setup
setup = Setup(name = "cylinders/4axons_flat")

# Geometry paramters
setup.geometry = Dict(
    :cell_shape => "cylinder",
    :ncell => 4,
    :rmin => 2.0,
    :rmax => 6.0,
    :dmin => 0.2,
    :dmax => 0.3,
    :deformation => (0.0, 0.0), #(0.05, 2π/3), #
    :height => 1.0,
    :include_in => true,
    :in_ratio => 0.6,
    :ecs_shape => "convex_hull",
    :ecs_ratio => 0.3,
    :refinement => 0.2,
)

# PDE paramaters
setup.pde = Dict(
    :σ_in => 0.002,
    :σ_out => 0.002,
    :σ_ecs => 0.002,
    :T₂_in => Inf,
    :T₂_out => Inf,
    :T₂_ecs => Inf,
    :ρ_in => 1.0,
    :ρ_out => 1.0,
    :ρ_ecs => 1.0,
    :κ_in_out => 1e-4,
    :κ_out_ecs => 1e-4,
    :κ_in => 0.0,
    :κ_out => 0.0,
    :κ_ecs => 0.0,
)

# Gradient sequences
setup.gradient = Dict(
    :directions => create_directions([1.0; 0.0; 0.0]),
    # :directions => create_directions(10, flat=true),
    :sequences => [PGSE(1000.0, 5000.0)],
    # :values => 0:500:2000,
    :values => [2000.0],
    :values_type => 'b',
)

# BTPDE solver
setup.btpde = Dict(:odesolver => ABDF2(), :reltol => 1e-4, :abstol => 1e-6, :nsave => 1)

# MF solver
setup.mf = Dict(:length_scale => 3, :neig_max => 400, :ninterval => 500)

# ImplicitEuler
# Trapezoid
# ABDF2
# QNDF
# QNDF1
# QNDF2
# QBDF
# QBDF1
# QBDF2
