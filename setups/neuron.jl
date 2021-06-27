## Neuron name
# name = "spindle/whole_neurons/03a_spindle2aFI"
name = "spindle/whole_neurons/03b_spindle4aACC"
# name = "spindle/whole_neurons/03b_spindle7aACC"
# name = "spindle/whole_neurons/04b_spindle3aFI"
# name = "spindle/whole_neurons/06b_spindle8aACC"
# name = "spindle/whole_neurons/16o_spindle13aFI"
# name = "spindle/whole_neurons/16o_spindle13aFI.stl"
# name = "spindle/separated_neurons/03a_spindle2aFI_soma"
# name = "spindle/separated_neurons/03a_spindle2aFI_dendrites_1"
# name = "spindle/separated_neurons/03a_spindle2aFI_dendrites_2"
# name = "spindle/separated_neurons/03b_spindle4aACC_dendrites_1"
# name = "spindle/separated_neurons/03b_spindle4aACC_dendrites_2"
# name = "spindle/separated_neurons/03b_spindle4aACC_soma"
# name = "spindle/separated_neurons/03b_spindle6aACC_dendrites_1"
# name = "spindle/separated_neurons/03b_spindle6aACC_dendrites_2"
# name = "spindle/separated_neurons/03b_spindle6aACC_soma"
# name = "spindle/separated_neurons/03b_spindle7aACC_dendrites_1"
# name = "spindle/separated_neurons/03b_spindle7aACC_dendrites_2"
# name = "spindle/separated_neurons/03b_spindle7aACC_soma"
# name = "spindle/separated_neurons/04b_spindle3aFI_dendrites_1"
# name = "spindle/separated_neurons/04b_spindle3aFI_dendrites_2"
# name = "spindle/separated_neurons/04b_spindle3aFI_soma"
# name = "spindle/separated_neurons/19o_spindle14aFI_dendrites_2.stl"
# name = "pyramidal/separated_neurons/02b_pyramidal1aACC_dendrites"
# name = "pyramidal/separated_neurons/02b_pyramidal1aACC_soma"
# name = "pyramidal/whole_neurons/02b_pyramidal1aACC"
# name = "pyramidal/whole_neurons/02a_pyramidal2aFI"
# name = "pyramidal/whole_neurons/04a_pyramidal4aACC"
# name = "pyramidal/whole_neurons/04a_pyramidal5aACC"
# name = "pyramidal/whole_neurons/04b_pyramidal6aACC"
# name = "pyramidal/whole_neurons/04b_pyramidal6aFI"
# name = "pyramidal/whole_neurons/25o_pyramidal18aFI"


## Setup
setup = NeuronSetup(
    name = name,
    ecs_shape = "no_ecs",
    ecs_ratio = 0.3,
    # refinement = 0.5,
    D_in = 0.002 * I(3),
    D_out = 0.002 * I(3),
    D_ecs = 0.002 * I(3),
    T₂_in = Inf,
    T₂_out = Inf,
    T₂_ecs = Inf,
    ρ_in = 1.0,
    ρ_out = 1.0,
    ρ_ecs = 1.0,
    κ_in_out = 1e-3,
    κ_out_ecs = 1e-3,
    κ_in = 0.0,
    κ_out = 0.0,
    κ_ecs = 0.0,
)

experiment = Experiment(
    gradient = (
        directions = create_directions([1.0; 0.0; 0.0]),
        sequences = [PGSE(2500.0, 4000.0)],
        values = [1000.0],
        values_type = "b",
    ),
    btpde = (odesolver = QNDF(), reltol = 1e-4, abstol = 1e-6, nsave = 1),
    btpde_midpoint = (θ = 0.5, timestep = 5),
    mf = (length_scale = 3, neig_max = 400, ninterval = 500),
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
# ROS3P(autodiff = false)
# Rosenbrock23(autodiff = false)
# ROS34PW1a(autodiff = false)
# ROS34PW3(autodiff = false)
# Rodas4(autodiff = false)
# Rodas4P(autodiff = false)
# Rodas5(autodiff = false)
# Cash4()

# Fixed time step
# MEBDF2()
