# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/SpinDoctor.jl")
    using .SpinDoctor
end

using SpinDoctor
using LinearAlgebra
using GLMakie

## Chose a plotting theme
set_theme!(theme_light())
set_theme!(theme_dark())
set_theme!(theme_black())


## Create model from setup recipe
# include("setups/axon.jl")
# include("setups/sphere.jl")
# include("setups/plates.jl")
# include("setups/cylinders.jl")
# include("setups/spheres.jl")
include("setups/neuron.jl")

mesh, = @time create_geometry(setup; recreate = true);
model = Model(; mesh, coeffs...);
volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)
ncompartment = length(model.mesh.points)
@info "Number of nodes per compartment:" length.(model.mesh.points)

## Plot mesh
plot_mesh(model.mesh)

## Assemble finite element matrices
matrices = @time assemble_matrices(model);


## Magnetic field gradient
TE = 5000.0
φ = -π / 6
R = [cos(φ) sin(φ) 0; -sin(φ) cos(φ) 0; 0 0 1]
g⃗(t) = 1.0 * R * [sin(2π * t / TE), sin(20π * t / TE) / 5, cos(2π * t / TE)]
gradient = GeneralGradient(; g⃗, TE)


## Solve BTPDE

# Callbacks for time stepping (plot solution, save time series)
printer = Printer(; nupdate = 1, verbosity = 2)
writer = VTKWriter(; nupdate = 5)
plotter = Plotter{T}(; nupdate = 1)
callbacks = [printer, plotter]
# callbacks = [printer, plotter, writer]

# General BTPDE for all gradients
btpde = BTPDE(; model, matrices)

# Solve BTPDE
ξ = @time solve(btpde, gradient; callbacks)


## Plot magnetization
plot_field(model.mesh, ξ)

## Compute signal
compute_signal(matrices.M, ξ)
compute_signal.(matrices.M_cmpts, split_field(model.mesh, ξ))


## Matrix Formalism

# Perform Laplace eigendecomposition
laplace = Laplace{T}(; model, matrices, neig_max = 400)
lap_eig = @time solve(laplace)
length_scales = eig2length.(lap_eig.values, D_avg)

# Truncate basis at minimum length scale
length_scale = 3
λ_max = length2eig(length_scale, D_avg)
lap_eig = limit_lengthscale(lap_eig, λ_max)

# Compute magnetization using the matrix formalism reduced order model
mf = MatrixFormalism(; model, matrices, lap_eig)
ξ = @time solve(mf, gradient; ninterval = 500)
