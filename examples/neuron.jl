# # [Neurons](@id neuron)
#
# The [RealNeuronMeshes](https://github.com/SpinDoctorMRI/RealNeuronMeshes) collection
# contains finite element meshes of 36 pyramidal neurons and 29 spindle neurons, created
# from the morphologies of [NeuroMorpho.org](https://neuromorpho.org) [Fang2020](@cite).
# Here, we compute the signal and the apparent diffusion coefficients of a spindle neuron.

using SpinDoctor
using LinearAlgebra
using CairoMakie
using Downloads
using p7zip_jll

# ## Mesh
#
# The meshes are zipped files in the Gmsh format. We download the spindle neuron
# `03b_spindle4aACC` and unzip it.

name = "03b_spindle4aACC"
url = "https://github.com/SpinDoctorMRI/RealNeuronMeshes/raw/master/volume_meshes/spindles/$name.msh.zip"
dir = mktempdir()
zipfile = Downloads.download(url, joinpath(dir, "$name.msh.zip"))
run(pipeline(`$(p7zip()) x $zipfile -o$dir -y`; stdout = devnull))
meshfile = joinpath(dir, "$name.msh")

# The [`NeuronSetup`](@ref) reads the mesh. The neuron could also be put in an
# extracellular space, as for the other recipes.

setup = NeuronSetup(; meshfile)
mesh, = create_geometry(setup)
plot_mesh(mesh; colors = [:mediumpurple], decorations = false, size = (700, 500))

# ## Signal
#
# We assume the same diffusion coefficient as in free water, and an impermeable membrane.

coeffs = coefficients(
    setup;
    D = (; neuron = 0.002 * I(3), ecs = 0.002 * I(3)),
    T2 = (; neuron = Inf, ecs = Inf),
    rho = (; neuron = 1.0, ecs = 1.0),
    kappa = (; neuron_ecs = 0.0, neuron = 0.0, ecs = 0.0),
    gamma = 2.67513e-4,
)
model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)
volume = only(matrices.volumes)

# The signal for gradients in the three axis directions, for b-values up to 4000 s/mm²:

profile = PGSE(10000.0, 40000.0)
bvalues = 0:500:4000
directions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
solver = ThetaSolver(; timestep = 200.0)
signals = map(directions) do dir
    map(bvalues) do b
        gradient = ScalarGradient(dir, profile; bvalue = b, gamma = model.gamma)
        u = solve_btpde(model, matrices, gradient, solver)
        abs(compute_signal(matrices.M, u)) / volume
    end
end

fig = Figure(; size = (600, 400))
ax = Axis(fig[1, 1]; xlabel = "b [s/mm²]", ylabel = "S / S₀", yscale = log10)
for (signal, label) in zip(signals, ["x", "y", "z"])
    scatterlines!(ax, bvalues, signal; label = "Gradient along $label")
end
lines!(
    ax,
    bvalues,
    exp.(-0.002 * bvalues);
    color = :gray,
    linestyle = :dash,
    label = "Free diffusion",
)
axislegend(ax; position = :lb)
fig

# The neuron restricts the diffusion, most across its dendrites. The magnetization for
# b = 4000 s/mm² along `x` shows where the spins dephase:

gradient = ScalarGradient(directions[1], profile; bvalue = 4000, gamma = model.gamma)
u = solve_btpde(model, matrices, gradient, solver)
plot_field(mesh, u; decorations = false, size = (750, 500))

# ## Apparent diffusion coefficients
#
# The homogenized ADC model gives the apparent diffusion tensor for all directions at once.
# Its eigenvalues show the anisotropy of the neuron:

tensor = only(solve_hadc_tensors(model, matrices, profile))
eigvals(Symmetric(tensor))

# The ADC in many directions, as a glyph:

directions = unitsphere(400)
adcs = [d' * tensor * d for d in eachcol(directions)]
plot_hardi(directions, adcs; label = "ADC [µm²/µs]")
