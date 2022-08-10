# The neuron meshes are provided in a separate git repository. Make sure to clone
#
#     https://github.com/SpinDoctorMRI/RealNeuronMeshes.git
#
# in order to run SpinDoctor simulations on real neuron geometries. You may need to manually
# unzip the mesh files of interest, as they are stored in text format.
#
# The file name is given by `prefix * name * postfix` (uncomment the desired file parts).

neuron_dir = "RealNeuronMeshes/volume_meshes"
isdir(neuron_dir) || error("""Cannot find neuron meshes.
                           Run `git clone https://github.com/SpinDoctorMRI/RealNeuronMeshes.git`
                           to download the mesh files.""")

## Choose whether to use whole neuron or certain parts
# prefix, postfix = "separated_", "_soma"
# prefix, postfix = "separated_", "_dendrites"
prefix, postfix = "separated_", "_dendrites_1"
# prefix, postfix = "separated_", "_dendrites_2"
# prefix, postfix = "", ""

## Neuron name
# name = "spindles/03a_spindle2aFI"
# name = "spindles/03a_spindle6aFI"
# name = "spindles/03b_spindle4aACC"
# name = "spindles/03b_spindle5aACC"
# name = "spindles/03b_spindle6aACC"
# name = "spindles/03b_spindle7aACC"
# name = "spindles/04b_spindle3aFI"
# name = "spindles/05b_spindle5aFI"
name = "spindles/06b_spindle8aACC"
# name = "spindles/07b_spindle9aACC"
# name = "spindles/08a_spindle13aACC"
# name = "spindles/09o_spindle7aFI"
# name = "spindles/09o_spindle8aFI"
# name = "spindles/10a_spindle18aACC"
# name = "spindles/12a_spindle19aACC"
# name = "spindles/12o_spindle9aFI"
# name = "spindles/13o_spindle10aFI"
# name = "spindles/15o_spindle12aFI"
# name = "spindles/16o_spindle13aFI"
# name = "spindles/19o_spindle14aFI"
# name = "spindles/21o_spindle15aFI"
# name = "spindles/23o_spindle16aFI"
# name = "spindles/25o_spindle17aFI"
# name = "spindles/26o_spindle18aFI"
# name = "spindles/27o_spindle19aFI"
# name = "spindles/28o_spindle20aFI"
# name = "spindles/28o_spindle21aFI"
# name = "spindles/29o_spindle22aFI"
# name = "spindles/30o_spindle23aFI"
# name = "pyramidals/02a_pyramidal2aFI"
# name = "pyramidals/02b_pyramidal1aACC"
# name = "pyramidals/02b_pyramidal1aFI"
# name = "pyramidals/03a_pyramidal9aFI"
# name = "pyramidals/03b_pyramidal2aACC"
# name = "pyramidals/03b_pyramidal3aACC"
# name = "pyramidals/03b_pyramidal3aFI"
# name = "pyramidals/03b_pyramidal4aFI"
# name = "pyramidals/03b_pyramidal9aFI"
# name = "pyramidals/04a_pyramidal4aACC"
# name = "pyramidals/04a_pyramidal5aACC"
# name = "pyramidals/04b_pyramidal5aFI"
# name = "pyramidals/04b_pyramidal6aACC"
# name = "pyramidals/04b_pyramidal6aFI"
# name = "pyramidals/04b_pyramidal7aACC"
# name = "pyramidals/05a_pyramidal10aACC"
# name = "pyramidals/05a_pyramidal8aACC"
# name = "pyramidals/05b_pyramidal7aFI"
# name = "pyramidals/05b_pyramidal8aFI"
# name = "pyramidals/05b_pyramidal9aACC"
# name = "pyramidals/06a_pyramidal11aACC"
# name = "pyramidals/06b_pyramidal10aFI"
# name = "pyramidals/06b_pyramidal12aACC"
# name = "pyramidals/07a_pyramidal13aACC"
# name = "pyramidals/07b_pyramidal14aACC"
# name = "pyramidals/08o_pyramidal11aFI"
# name = "pyramidals/10a_pyramidal15aACC"
# name = "pyramidals/11a_pyramidal16aACC"
# name = "pyramidals/11o_pyramidal12aFI"
# name = "pyramidals/17o_pyramidal13aFI"
# name = "pyramidals/18o_pyramidal14aFI"
# name = "pyramidals/20o_pyramidal15aFI"
# name = "pyramidals/22o_pyramidal16aFI"
# name = "pyramidals/24o_pyramidal17aFI"
# name = "pyramidals/25o_pyramidal18aFI"
# name = "pyramidals/31o_pyramidal19aFI"

## Resulting neuron filename
name = prefix * name * postfix
meshdir = joinpath(neuron_dir, name)

# Floating point type for simulations
T = Float64

# Geometrical setup
setup = NeuronSetup{T}(;
    ecs = NoECS(),
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
