# # [Dense tissue with periodic boundary conditions](@id periodic)
#
# A simulated sample of tissue has boundaries that the tissue does not have. With
# impermeable faces, the water of the extracellular space (ECS) is trapped in the sample,
# and its apparent diffusion coefficient (ADC) decays to zero at long diffusion times. In
# the tissue, the ECS is connected, and its ADC tends to the free diffusion coefficient
# divided by the squared tortuosity. Periodic boundary conditions make the sample a
# periodic tissue, without walls.

using SpinDoctor
using LinearAlgebra
using Random
using CairoMakie

# ## Sample
#
# Ten polyhedral cells ([`VoronoiSetup`](@ref)) fill 80 % of a periodic box of 12 µm. The
# cells are cut by the faces of the box, and their periodic images enter on the opposite
# faces. The mesh is periodic: the opposite faces have the same points.

ncell = 10
setup = VoronoiSetup(; ncell, boxsize = 12.0, fraction = 0.8)
mesh, surfaces, cells = create_geometry(setup; rng = Xoshiro(1))
plot_mesh(mesh; colors = [fill(:seagreen, ncell); :lightgray], decorations = false)

# The ECS alone is a network of thin sheets between the cells:

plot_mesh(mesh; compartments = [ncell + 1], colors = fill(:lightsteelblue, ncell + 1))

# The cells are impermeable, with the same diffusion coefficient as the ECS. The
# coefficients are given as for spheres with one layer.

coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3)], ecs = 0.002 * I(3)),
    T2 = (; cell = [Inf], ecs = Inf),
    rho = (; cell = [1.0], ecs = 1.0),
    kappa = (; cell_interfaces = Float64[], cell_boundaries = [0], cell_ecs = 0, ecs = 0),
    gamma = 2.67513e-4,
)
model = Model(; mesh, coeffs...);

# The matrices with impermeable box faces, and with periodic boundary conditions:

walls = assemble_matrices(model)
periodic = assemble_matrices(model; periodic = true);

# ## ADC and diffusion time
#
# The homogenized ADC model gives the apparent diffusion tensors of the compartments for
# PGSE sequences with short pulses (1 ms), for diffusion times between 2 and 100 ms.

Deltas = [2000.0, 5000.0, 10000.0, 20000.0, 50000.0, 100000.0]
D0 = 0.002
mean_adc(tensor) = tr(tensor) / 3
adcs = map((; walls, periodic)) do matrices
    map(Deltas) do Delta
        profile = PGSE(1000.0, Delta)
        tensors = solve_hadc_tensors(model, matrices, profile; timestep = Delta / 200)
        mean_adc.(tensors) ./ D0
    end
end;

# With impermeable faces, the ADC of the ECS decays at long diffusion times, as for a
# closed pore. With periodic boundary conditions, the ECS is a connected network of sheets,
# much thinner than the distance the water diffuses, and its ADC is already at its long-time
# limit D₀ / τ², which defines the tortuosity τ of the ECS. Sheets with random orientations
# would give 2/3, the junctions between them increase it. The ADC of the closed cells
# decays in both cases (the cells cut by the box faces differ, as periodicity joins their
# pieces).

fig = Figure(; size = (650, 420))
ax = Axis(
    fig[1, 1];
    xlabel = "Diffusion time Δ [ms]",
    ylabel = "ADC / D₀",
    xscale = log10,
    limits = (nothing, (0, 1)),
)
Δ = Deltas ./ 1000
scatterlines!(ax, Δ, [a[end] for a in adcs.periodic]; label = "ECS, periodic")
scatterlines!(ax, Δ, [a[end] for a in adcs.walls]; label = "ECS, impermeable faces")
scatterlines!(
    ax,
    Δ,
    [sum(a[1:ncell]) / ncell for a in adcs.periodic];
    label = "Cells (mean)",
)
axislegend(ax; position = :rt)
fig

# The tortuosity of the ECS, from the longest diffusion time:

tortuosity = 1 / sqrt(adcs.periodic[end][end])
