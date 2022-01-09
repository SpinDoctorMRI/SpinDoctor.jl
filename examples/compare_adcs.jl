using LinearAlgebra
using GLMakie

# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/SpinDoctor.jl")
    using .SpinDoctor
else
    using SpinDoctor
end

set_theme!(theme_dark())


## Choose setup recipe
# include("setups/axon.jl")
# include("setups/sphere.jl")
include("setups/plates.jl")
# include("setups/cylinders.jl")
# include("setups/spheres.jl")
# include("setups/neuron.jl")

## Create model from setup recipe
femesh, = @time create_geometry(setup; recreate = true)
model = Model(; mesh = femesh, coeffs...)
volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)
ncompartment = length(model.mesh.points)
matrices = @time assemble_matrices(model);
@info "Number of nodes per compartment:" length.(model.mesh.points)


## Magnetic field gradient
dir = [1.0, 1.0, 1.0]
profile = PGSE(2500.0, 4000.0)
b = 1000
g = √(b / int_F²(profile)) / model.γ
gradient = ScalarGradient(dir, profile, g)

## Short term approximation (STA)
adc_sta_cmpts = compute_adc_sta(model, gradient)

## Solve BTPDE and fit ADC for different b-values
bvalues = 0:100:1000
gvalues = map(b -> √(b / int_F²(profile)) / coeffs.γ, bvalues)
gradients = [ScalarGradient(gradient.dir, gradient.profile, g) for g ∈ gvalues]
# btpde = GeneralBTPDE(;
#     model, matrices, reltol = 1e-4, abstol = 1e-6,
# )
btpde = IntervalConstanBTPDE{T}(; model, matrices, θ = 0.5, timestep = 5)
ξ, = solve_multigrad(btpde, gradients)
signals = [compute_signal(matrices.M, ξ) for ξ ∈ ξ]
signals_cmpts = [compute_signal.(matrices.M_cmpts, split_field(model.mesh, ξ)) for ξ ∈ ξ]
adc_fit = fit_adc(bvalues, signals)
adc_fit_cmpts =
    [fit_adc(bvalues, [s[icmpt] for s ∈ signals_cmpts]) for icmpt = 1:ncompartment]


## Solve HADC
hadc = HADC(; model, matrices, reltol = 1e-4, abstol = 1e-6)
adc_homogenized_cmpts = @time solve(hadc, gradient)


## Matrix Formalism ADC
# Perform Laplace eigendecomposition
laplace = Laplace{T}(; model, matrices, neig_max = 400)
lap_eig = @time solve(laplace)

D_mf = compute_mf_diffusion_tensor(model.mesh, matrices.M, lap_eig, gradient)
adc_mf = dir'D_mf * dir / dir'dir
