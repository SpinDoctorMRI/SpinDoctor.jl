# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/SpinDoctor.jl")
    using .SpinDoctor
end

using SpinDoctor
using LinearAlgebra
using GLMakie

## Create model from setup recipe
# include("setups/axon.jl")
# include("setups/cylinders.jl")
# include("setups/disks.jl")
# include("setups/neuron.jl")
include("setups/plates.jl")
# include("setups/slabs.jl")
# include("setups/sphere.jl")
# include("setups/spheres.jl")

mesh, surfaces, cells = create_geometry(setup; savedir, recreate = true);
model = Model(; mesh, coeffs...);
dim = size(surfaces.points, 1)
@info "Number of nodes per compartment:" length.(model.mesh.points)

## Plot mesh
plot_surfaces(surfaces)
plot_surfaces(surfaces, 1:10)
plot_mesh(mesh, 1:4, 10:10)
plot_mesh(mesh)

## Assemble finite element matrices
matrices = assemble_matrices(model);

## Magnetic field gradient
if dim == 2
    dir = [1.0, 0.0]
elseif dim == 3
    dir = [1.0, 0.0, 0.0]
end
profile = PGSE(2000.0, 6000.0)
# profile = CosOGSE(5000.0, 5000.0, 2)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
gradient = ScalarGradient(dir, profile, g)

## Solve BTPDE

# Callbacks for time stepping (plot solution, save time series)
printer = Printer(; nupdate = 1, verbosity = 2)
writer = VTKWriter(; dir = joinpath("output", name), nupdate = 5)
plotter = Plotter(; nupdate = 5)
callbacks = [printer, plotter]
# callbacks = [printer, plotter, writer]

# Choose BTDPE solver (specialized solver only for PGSE)
solver = IntervalConstantSolver(; θ = 0.5, timestep = 5.0)
solver = QNDF()

# Solve BTPDE
btpde = BTPDE(; model, matrices)
ξ = @time solve(btpde, gradient, solver; callbacks)

## Plot magnetization
plot_field(model.mesh, ξ)

## Compute signal
compute_signal(matrices.M, ξ)
compute_signal.(matrices.M_cmpts, split_field(model.mesh, ξ))

## Save magnetization
savefield(model.mesh, ξ, joinpath("output", name, "magnetization"))

## Matrix Formalism

volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / dim * tr.(model.D)' * volumes / sum(volumes)

# Perform Laplace eigendecomposition
laplace = Laplace(; model, matrices, neig_max = 400)
lap_eig = @time solve(laplace)
length_scales = eig2length.(lap_eig.values, D_avg)

# Truncate basis at minimum length scale
length_scale = 3
λ_max = length2eig(length_scale, D_avg)
lap_eig = limit_lengthscale(lap_eig, λ_max)

# Compute magnetization using the matrix formalism reduced order model
mf = MatrixFormalism(; model, matrices, lap_eig)
ξ = @time solve(mf, gradient; ninterval = 500)

## Solve HADC
hadc = HADC(; model, matrices)
adc_cmpts = @time solve(hadc, gradient)

## Solve Karger model

# Compute HADC and fit difftensors
directions = unitsphere(50)
gradients = [
    ScalarGradient(collect(d), gradient.profile, gradient.amplitude) for
    d ∈ eachcol(directions)
]
hadc = HADC(; model, matrices)
adcs = Vector{eltype(directions)}[]
@time for grad ∈ gradients
    @show grad.dir
    push!(adcs, solve(hadc, grad))
end
difftensors = fit_tensors(directions, adcs)

# Solve Karger
karger = Karger(; model, difftensors)
signal = @time solve(karger, gradient; timestep = 5.0)

## Solve analytical model
# Compute analytical Laplace eigenfunctions
length_scale = 0.3
eigstep = 1e-8
eiglim = length2eig(length_scale, D_avg)
analytical_coeffs = analytical_coefficients(setup, coeffs)
analytical_laplace = AnalyticalLaplace(; analytical_coeffs..., eiglim, eigstep)
lap_mat = @time solve(analytical_laplace)

# Compute analytical matrix formalism signal truncation
volumes = get_cmpt_volumes(model.mesh)
analytical_mf = AnalyticalMatrixFormalism(; analytical_laplace, lap_mat, volumes)
signal = solve(analytical_mf, gradient)
