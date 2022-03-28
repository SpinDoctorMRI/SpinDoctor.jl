# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/SpinDoctor.jl")
    using .SpinDoctor
end
using MKL
using SpinDoctor
using LinearAlgebra
using CairoMakie

## Chose a plotting theme
set_theme!(theme_light())
set_theme!(theme_dark())
set_theme!(theme_black())


## Create model from setup recipe
# include("setups/axon.jl")
# include("setups/sphere.jl")
# include("setups/plates.jl")
include("setups/cylinders.jl")
# include("setups/spheres.jl")
# include("setups/neuron.jl")

model, matrices, surfaces, = prepare_simulation(setup;recreate = false)

savepath = add_checksum("res_mf",setup)

## Plot mesh
plot_surfaces(surfaces, 1:3)
plot_mesh(model.mesh, 1:1)

## Magnetic field gradient
dir = [1.0, 0.0, 0.0]
profile = PGSE(2000.0, 6000.0)
# profile = CosOGSE(5000.0, 5000.0, 2)
b = 1000
g = √(b / int_F²(profile)) / model.γ
gradient = ScalarGradient(dir, profile, g)


# ## Solve BTPDE

# # Callbacks for time stepping (plot solution, save time series)
# printer = Printer(; nupdate = 1, verbosity = 2)
# writer = VTKWriter(; nupdate = 5)
# plotter = Plotter{T}(; nupdate = 5)
# # callbacks = [printer, plotter]
# callbacks = [printer, plotter, writer]

# # Choose BTDPE solver (specialized solver only for PGSE)
# solver = IntervalConstantSolver{T}(; θ = 0.5, timestep = 5.0)
# solver = QNDF(autodiff = false)

# # Solve BTPDE
# btpde = BTPDE(; model, matrices)
# ξ = @time solve(btpde, gradient, solver; callbacks)


# ## Plot magnetization
# plot_field(model.mesh, ξ)

# ## Compute signal
# S = compute_signal(matrices.M, ξ)
# S_cmpts = compute_signal.(matrices.M_cmpts, split_field(model.mesh, ξ))
# @info S, S_cmpts, ξ
# ## Save magnetization
# savefield(model.mesh, ξ, "output/magnetization")


## Matrix Formalism

# Perform Laplace eigendecomposition
laplace = Laplace{T}(; model, matrices, neig_max = 400, length_scale = 3)
lap_eig = @time solve(laplace;rerun=true,savepath=savepath)

# Compute magnetization using the matrix formalism reduced order model
mf = MatrixFormalism(; model, matrices, lap_eig)
S, S_cmpts, ξ = @time solve(mf, gradient; ninterval = 500)
lap_eig.λ[end], lap_eig.ϕ[end]


# ## Solve HADC
# hadc = HADC(; model, matrices)
# adc_cmpts = @time solve(hadc, gradient)


## Solve Karger model

# # Compute HADC and fit difftensors
# directions = unitsphere(50)
# gradients = [
#     ScalarGradient(collect(d), gradient.profile, gradient.amplitude) for
#     d ∈ eachcol(directions)
# ]
# hadc = HADC(; model, matrices)
# adcs, = @time solve_multigrad(hadc, gradients)
# difftensors = fit_tensors(directions, adcs)

# # Solve Karger
# karger = Karger(; model, difftensors)
# signal = @time solve(karger, gradient; timestep = 5.0)


# ## Solve analytical model
# # Compute analytical Laplace eigenfunctions
# length_scale = 0.3
# eigstep = 1e-8
# eiglim = length2eig(length_scale, D_avg)
# analytical_coeffs = analytical_coefficients(setup, coeffs)
# analytical_laplace = AnalyticalLaplace(; analytical_coeffs..., eiglim, eigstep)
# lap_mat = @time solve(analytical_laplace)

# # Compute analytical matrix formalism signal truncation
# analytical_mf = AnalyticalMatrixFormalism(; analytical_laplace, lap_mat, volumes)
# signal = solve(analytical_mf, gradient)
