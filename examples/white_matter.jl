# # [Anisotropic diffusion in white matter](@id white-matter)
#
# In white matter, myelinated axons are packed densely along a common direction. Water
# diffuses freely along the axons, and is restricted across them. This example simulates a
# sample of such tissue in a periodic box, with periodic boundary conditions, so that the
# sample represents a larger piece of tissue.

using SpinDoctor
using LinearAlgebra
using Random
using CairoMakie

# ## Sample
#
# Fifteen myelinated axons (radii between 1 and 2.5 µm) fill 65 % of the cross-section
# ([`PeriodicBoxECS`](@ref)). The axons are cut by the faces of the box, and their
# periodic images enter on the opposite faces.

ncell = 15
setup = CylinderSetup(;
    ncell,
    rmin = 1.0,
    rmax = 2.5,
    nsidewall = 16,
    dmin = 0.15,
    layersizes = [0.7, 1.0],
    ecs = PeriodicBoxECS(; fraction = 0.65),
    height = 8.0,
)
mesh, surfaces, cells = create_geometry(setup; rng = Xoshiro(4))
colors = [fill(:steelblue, ncell); fill(:orange, ncell); :lightgray]
plot_mesh(mesh; colors, decorations = false)

# The axons, the myelin and the ECS exchange water through permeable membranes. Water
# diffuses more slowly in the myelin.

coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3), 0.0005 * I(3)], ecs = 0.002 * I(3)),
    T2 = (; cell = [Inf, Inf], ecs = Inf),
    rho = (; cell = [1.0, 1.0], ecs = 1.0),
    kappa = (;
        cell_interfaces = [1e-5],
        cell_boundaries = [0, 0],
        cell_ecs = 1e-5,
        ecs = 0,
    ),
    gamma = 2.67513e-4,
)
model = Model(; mesh, coeffs...);

# With `periodic = true`, the points on opposite faces of the box are identified, and the
# solvers use periodic boundary conditions instead of impermeable box faces.

matrices = assemble_matrices(model; periodic = true);

# ## Apparent diffusion tensor
#
# The homogenized ADC model gives the apparent diffusion tensor of each compartment for all
# directions at once, here for a PGSE sequence with a diffusion time of 20 ms. The tensor of
# the sample is the average of the compartments, weighted by their contributions to the
# signal ([`compartment_weights`](@ref)).

profile = PGSE(5000.0, 20000.0)
tensors = solve_hadc_tensors(model, matrices, profile)
weights = compartment_weights(model)
tensor = sum(weights .* tensors)

# The ADC is largest along the axons (`z`), with the free diffusion coefficients (the
# periodic boundary conditions do not restrict the diffusion along the axons), and smaller
# across:

directions = unitsphere(500)
adcs = [d' * tensor * d for d in eachcol(directions)]
plot_hardi(directions, adcs; label = "ADC [µm²/µs]", decorations = false)

# ## Signal attenuation
#
# The Bloch-Torrey equation gives the signal beyond the ADC regime. We compare gradients
# across (`x`) and along (`z`) the axons.

bvalues = [0, 500, 1000, 2000, 3000, 4000]
solver = ThetaSolver(; timestep = 200.0)
signals = map((; across = [1.0, 0.0, 0.0], along = [0.0, 0.0, 1.0])) do dir
    map(bvalues) do b
        gradient = ScalarGradient(dir, profile; bvalue = b, gamma = model.gamma)
        u = solve_btpde(model, matrices, gradient, solver)
        abs(compute_signal(matrices.M, u)) / sum(matrices.volumes)
    end
end

# Along the axons, each compartment would attenuate as in free diffusion (dashed line)
# without exchange. The exchange of water between the slow myelin and the fast compartments
# makes the slow component decay faster. Across the axons, the signal attenuates slowly,
# mostly from the ECS, and deviates from a monoexponential decay at higher b-values.

free = [sum(weights .* exp.(-b .* [D[3, 3] for D in model.D])) for b in bvalues]
fig = Figure(; size = (600, 400))
ax = Axis(fig[1, 1]; xlabel = "b [s/mm²]", ylabel = "S / S₀", yscale = log10)
scatterlines!(ax, bvalues, signals.across; label = "Across the axons")
scatterlines!(ax, bvalues, signals.along; label = "Along the axons")
lines!(ax, bvalues, free; color = :black, linestyle = :dash, label = "Free, no exchange")
axislegend(ax; position = :lb)
fig

# The magnetization across the axons, for b = 4000 s/mm²: the spins dephase in the ECS,
# while the water in the axons keeps its phase.

gradient = ScalarGradient([1.0, 0.0, 0.0], profile; bvalue = 4000, gamma = model.gamma)
u = solve_btpde(model, matrices, gradient, solver)
plot_field(mesh, u; decorations = false)

# The periodic time stepping runs much faster on an NVIDIA GPU, with
# `ThetaSolver(; timestep, arraytype = CuArray)` after `using CUDA, CUDSS`.
