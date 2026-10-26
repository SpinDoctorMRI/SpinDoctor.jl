# # [Apparent diffusion coefficients](@id adc-methods)
#
# At low b-values, the signal attenuation is described by the apparent diffusion
# coefficient (ADC), ``S / S_0 \approx e^{-b \, \mathrm{ADC}}`` (see
# [Apparent diffusion coefficient](@ref adc)). SpinDoctor computes the ADC in several ways:
#
# - by fitting the signals of the Bloch-Torrey equation for several b-values (the
#   reference, but the most expensive);
# - with the homogenized ADC model ([`solve_hadc`](@ref)), which assumes impermeable
#   membranes;
# - with the matrix formalism ([`compute_mf_diffusion_tensor`](@ref)), from the Laplace
#   eigenfunctions;
# - with the short diffusion time approximation ([`compute_adc_sta`](@ref)), from the
#   surface-to-volume ratios.
#
# In this example, we compare them for stacked slabs separated by permeable membranes.

using SpinDoctor
using LinearAlgebra
using CairoMakie

# ## Model
#
# Five slabs of 5 µm, twisted around the vertical axis. The water diffuses freely along
# the slabs, and crosses the membranes between them slowly. The second slab has a higher
# diffusion coefficient.

ncell = 5
setup = SlabSetup(;
    depth = 50.0,
    widths = fill(5.0, ncell),
    height = 50.0,
    twist = pi / 6,
    refinement = 10.0,
)
D = fill(0.002, ncell)
D[2] = 0.003
coeffs = coefficients(
    setup;
    D = [d * I(3) for d in D],
    T2 = fill(Inf, ncell),
    rho = fill(1.0, ncell),
    kappa = (; interfaces = fill(1e-4, ncell - 1), boundaries = fill(0.0, ncell)),
    gamma = 2.67513e-4,
)
mesh, = create_geometry(setup)
plot_mesh(mesh; decorations = false)

#-

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)
volumes = matrices.volumes
D_mean = sum(volumes .* D) / sum(volumes)

# The gradient is a PGSE sequence across and along the slabs (the `x` and `z` directions).

dir = [1.0, 0.0, 1.0]
profile = PGSE(2500.0, 4000.0)
gradient = ScalarGradient(dir, profile; bvalue = 1000, gamma = model.gamma)

# ## ADCs
#
# The signals of the Bloch-Torrey equation for b-values up to 4000 s/mm², and their ADC:
# the initial slope of a quadratic fit of the log-signal (for the sample, and for each
# compartment):

bvalues = 0:400:4000
solver = ThetaSolver(; timestep = 20.0)
signals = map(bvalues) do b
    gradient = ScalarGradient(dir, profile; bvalue = b, gamma = model.gamma)
    u = solve_btpde(model, matrices, gradient, solver)
    compute_signal.(matrices.M_cmpts, split_field(mesh, u))
end
adc_fit = fit_adc(bvalues, sum.(signals); degree = 2)
adc_fit_cmpts = [fit_adc(bvalues, [s[i] for s in signals]; degree = 2) for i = 1:ncell]

# The homogenized ADC model, with impermeable membranes:

adc_hadc_cmpts = solve_hadc(model, matrices, gradient)
adc_hadc = sum(volumes .* adc_hadc_cmpts) / sum(volumes)

# The matrix formalism, with the Laplace eigenfunctions of the whole sample (including the
# permeable membranes), gives the effective diffusion tensor:

laplace = solve_laplace(model, matrices; neig_max = 400)
D_mf = compute_mf_diffusion_tensor(mesh, matrices.M, laplace, gradient)
adc_mf = dir' * D_mf * dir / (dir' * dir)

# The short diffusion time approximation:

adc_sta_cmpts = compute_adc_sta(model, gradient)
adc_sta = sum(volumes .* adc_sta_cmpts) / sum(volumes)

# ## Comparison
#
# The ADCs of the sample, relative to the mean diffusion coefficient, and the signal
# attenuations they predict, with the signals of the Bloch-Torrey equation:

methods = ["Fit", "HADC", "MF", "STA", "Free"]
adcs = [adc_fit, adc_hadc, adc_mf, adc_sta, D_mean]
fig = Figure(; size = (900, 400))
ax = Axis(fig[1, 1]; xticks = (1:5, methods), ylabel = "ADC / D", title = "ADC")
barplot!(ax, 1:5, adcs ./ D_mean; color = Makie.wong_colors()[1:5])
ax = Axis(
    fig[1, 2];
    xlabel = "b [s/mm²]",
    ylabel = "S / S₀",
    yscale = log10,
    title = "Signal attenuation",
)
for (method, adc, color) in zip(methods, adcs, Makie.wong_colors())
    lines!(ax, [0, 4000], [1, exp(-4000 * adc)]; color, linestyle = :dash, label = method)
end
scatter!(ax, bvalues, abs.(sum.(signals)) ./ sum(volumes); color = :black, label = "BTPDE")
axislegend(ax; position = :lb)
fig

# - The signal deviates from a straight line (in log scale) at the higher b-values: the
#   ADC only describes the low b-values.
# - The HADC is smaller than the fitted ADC, as it assumes impermeable membranes. The MF
#   effective tensor, from the same model with permeable membranes, is close to the HADC
#   here, as the membranes are only weakly permeable.
# - The short time approximation underestimates the ADC, as the diffusion time is not
#   short.
# - Free diffusion attenuates much more: the membranes restrict the diffusion.

# The compartment ADCs: the fitted ADCs are larger than the HADC, as the water crosses the
# membranes, most for the inner slabs, which have membranes on both sides. The HADC and the
# STA treat each slab separately.

fig = Figure(; size = (650, 350))
ax = Axis(fig[1, 1]; xlabel = "Slab", ylabel = "ADC / D", xticks = 1:ncell)
group = repeat(1:3; inner = ncell)
values = [adc_fit_cmpts ./ D; adc_hadc_cmpts ./ D; adc_sta_cmpts ./ D]
barplot!(
    ax,
    repeat(1:ncell, 3),
    values;
    dodge = group,
    color = Makie.wong_colors()[[1, 2, 4]][group],
)
labels = ["Fit", "HADC", "STA"]
elements = [PolyElement(; color = Makie.wong_colors()[i]) for i in [1, 2, 4]]
Legend(fig[1, 2], elements, labels)
fig
