## Neuron name
# name = "spindle/whole_neurons/03a_spindle2aFI"
# name = "spindle/whole_neurons/03b_spindle4aACC"
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
name = "spindle/separated_neurons/03b_spindle6aACC_soma"
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

neuron_dir = "NeuronMeshes"
isdir(neuron_dir) || error("""Cannot find neuron meshes.
                           Run `git clone https://github.com/SpinDoctorMRI/NeuronMeshes.git`
                           to download the mesh files.""")

# Floating point type for simulations
T = Float64

# Geometrical setup
setup = NeuronSetup{T}(;
    name,
    meshdir = neuron_dir,
    ecs_shape = :no_ecs,
    ecs_ratio = 0.3,
    # refinement = 0.5,
)

# Get compartimentalized coefficient vectors
coeffs = coefficients(
    setup;
    D = (; neuron = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; neuron = Inf, ecs = Inf),
    ρ = (; neuron = 1.0, ecs = 1.0),
    κ = (; neuron_ecs = 1e-4, in = 0.0, neuron = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)
