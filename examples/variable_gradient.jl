using LinearAlgebra
using GLMakie
using OrdinaryDiffEq: Rodas4

# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/SpinDoctor.jl")
    using .SpinDoctor
else
    using SpinDoctor
end

set_theme!(theme_dark())


## Create model from setup recipe
# include("setups/axon.jl")
# include("setups/sphere.jl")
# include("setups/cylinders.jl")
# include("setups/spheres.jl")
include("setups/neuron.jl")

femesh, = @time create_geometry(setup)
model = Model(; mesh = femesh, coeffs...)
volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)
ncompartment = length(model.mesh.points)
@info "Number of nodes per compartment:" length.(model.mesh.points)

## Plot mesh
plot_mesh(model.mesh)

## Assemble finite element matrices
matrices = @time assemble_matrices(model);


## Magnetic field gradient
TE = 5000
g⃗(t) = 1.0 * [sin(2π * t / TE), sin(20π * t / TE) / 5, cos(2π * t / TE)]
gradient = GeneralGradient{T,typeof(g⃗)}(; g⃗, TE)



## Solve BTPDE

# Callbacks for time stepping (plot solution, save time series)
printer = Printer(; nupdate = 1, verbosity = 2)
writer = VTKWriter(; nupdate = 5)
plotter = Plotter{T}(; nupdate = 5)
callbacks = [printer, plotter]
# callbacks = [printer, plotter, writer]

# General BTPDE for all gradients
btpde = GeneralBTPDE(;
    model, matrices, reltol = 1e-4, abstol = 1e-6, odesolver = Rodas4(autodiff = false)
)

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
mf = MatrixFormalism(; model, matrices, lap_eig, ninterval = 500)
ξ = @time solve(mf, gradient)
