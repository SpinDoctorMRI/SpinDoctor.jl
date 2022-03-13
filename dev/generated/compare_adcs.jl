using SpinDoctor
using LinearAlgebra

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

ncell = 5
setup = PlateSetup(;
    name = "Plates",
    width = 50.0,
    depth = 50.0,
    heights = fill(5.0, ncell),
    bend = 0.0,
    twist = 0.0,
    refinement = 10.0,
)
coeffs = coefficients(
    setup;
    D = [0.002 * I(3) for _ = 1:ncell],
    T₂ = fill(Inf, ncell),
    ρ = fill(1.0, ncell),
    κ = (; interfaces = fill(1e-4, ncell - 1), boundaries = fill(0.0, ncell)),
    γ = 2.67513e-4,
)

mesh, = create_geometry(setup; recreate = true)
plot_mesh(mesh)

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model);

volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)
ncompartment = length(model.mesh.points)

dir = [1.0, 0.0, 1.0]
profile = PGSE(2500.0, 4000.0)
b = 1000
g = √(b / int_F²(profile)) / model.γ
gradient = ScalarGradient(dir, profile, g)

adc_sta_cmpts = compute_adc_sta(model, gradient)
adc_sta = volumes'adc_sta_cmpts / sum(volumes)

bvalues = 0:400:4000
gvalues = map(b -> √(b / int_F²(profile)) / coeffs.γ, bvalues)
gradients = [ScalarGradient(gradient.dir, gradient.profile, g) for g ∈ gvalues]

btpde = BTPDE(; model, matrices)
solver = IntervalConstantSolver(; θ = 0.5, timestep = 5.0)
ξ, = solve_multigrad(btpde, gradients, solver)

signals = [compute_signal(matrices.M, ξ) for ξ ∈ ξ]
signals_cmpts = [compute_signal.(matrices.M_cmpts, split_field(model.mesh, ξ)) for ξ ∈ ξ]

adc_fit = fit_adc(bvalues, signals)
adc_fit_cmpts =
    [fit_adc(bvalues, [s[icmpt] for s ∈ signals_cmpts]) for icmpt = 1:ncompartment]

hadc = HADC(; model, matrices)
adc_homogenized_cmpts = solve(hadc, gradient)
adc_homogenized = volumes'adc_homogenized_cmpts / sum(volumes)

laplace = Laplace(; model, matrices, neig_max = 400)
lap_eig = solve(laplace)

D_mf = compute_mf_diffusion_tensor(model.mesh, matrices.M, lap_eig, gradient)

adc_mf = dir'D_mf * dir / dir'dir

n = ncompartment
fig = Figure()
ax = Axis(fig[1, 1];
    xticks = (1:4, ["STA", "Fit BTPDE", "HADC", "MF"]),
    ylabel = "ADC / D",
    title = "Compartment ADCs",
)
barplot!(ax, fill(1, n), adc_sta_cmpts ./ D_avg; dodge = 1:n)
barplot!(ax, fill(2, n), adc_fit_cmpts ./ D_avg; dodge = 1:n)
barplot!(ax, fill(3, n), adc_homogenized_cmpts ./ D_avg; dodge = 1:n)
barplot!(ax, [4], [adc_mf / D_avg])
fig

fig = Figure()
ax = Axis(fig[1,1]; xlabel = "b", yscale = log10, title = "Signal attenuation")
lines!(ax, [0, bvalues[end]], [1, exp(-adc_sta * bvalues[end])]; linestyle = :dash, label = "ADC STA")
lines!(ax, [0, bvalues[end]], [1, exp(-adc_fit * bvalues[end])]; linestyle = :dash, label = "ADC Fit")
lines!(ax, [0, bvalues[end]], [1, exp(-adc_homogenized * bvalues[end])]; linestyle = :dash, label = "HADC")
lines!(ax, [0, bvalues[end]], [1, exp(-adc_mf * bvalues[end])]; linestyle = :dash, label = "ADC MF")
lines!(ax, [0, bvalues[end]], [1, exp(-D_avg * bvalues[end])]; linestyle = :dash, label = "Free diffusion")
scatterlines!(ax, bvalues, abs.(signals) ./ abs(signals[1]); label = "BTPDE Signal")
axislegend(ax)
fig

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

