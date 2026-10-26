# # [Getting started: myelinated axons](@id quickstart)
#
# This example walks through a complete simulation: the geometry of three myelinated
# axons in an extracellular space (ECS), the finite element model, and the diffusion MRI
# signal for gradients across and along the axons.

using SpinDoctor
using LinearAlgebra
using Random
using CairoMakie # Or GLMakie, for interactive 3D plots

# ## Geometry
#
# The [`CylinderSetup`](@ref) recipe creates parallel cylinders from a random
# configuration of disks. Each cell has two layers: the axon (the inner 60 % of the radius)
# and the myelin sheath. The ECS fills the convex hull of the cells, enlarged by a margin.
# Lengths are in micrometers.

setup = CylinderSetup(;
    ncell = 3,
    rmin = 2.0,
    rmax = 4.0,
    nsidewall = 20,
    layersizes = [0.6, 1.0],
    ecs = ConvexHullECS(; margin = 1.0),
    height = 20.0,
)
mesh, surfaces, cells = create_geometry(setup; rng = Xoshiro(1));

# The compartments are numbered with the inner layers (the axons) first, then the outer
# layers (the myelin), and the ECS last. We color them accordingly, and make the ECS
# transparent.

ncell = 3
colors = [fill(:steelblue, ncell); fill(:orange, ncell); :gray]
alpha = [fill(1.0, 2ncell); 0.25]
plot_mesh(mesh; colors, alpha, decorations = false)

# ## Model
#
# The coefficients are given for the layers of the cells and for the ECS: diffusion
# tensors `D` (µm²/µs), T2 relaxation times `T2` (µs), initial spin densities `rho`, and
# permeabilities `kappa` (µm/µs) of the interfaces between the layers, of the outer
# boundaries of the cells, and between the cells and the ECS. Water diffuses more slowly in
# the myelin.

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
model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model);

# ## Diffusion MRI signal
#
# A pulsed gradient spin echo (PGSE) sequence has two pulses of duration `delta` separated
# by `Delta` (µs). We solve the Bloch-Torrey equation for the magnetization at the echo
# time, for gradients across the axons (`x`) and along them (`z`), with b = 1000 s/mm².

profile = PGSE(5000.0, 20000.0)
solver = ThetaSolver(; timestep = 50.0)
directions = (; across = [1.0, 0.0, 0.0], along = [0.0, 0.0, 1.0])
magnetizations = map(directions) do dir
    gradient = ScalarGradient(dir, profile; bvalue = 1000, gamma = model.gamma)
    solve_btpde(model, matrices, gradient, solver)
end;

# The signal is the integral of the magnetization, here for each compartment, relative to
# the signal without gradient (the volume of the compartment).

attenuation = map(magnetizations) do u
    signals = compute_signal.(matrices.M_cmpts, split_field(mesh, u))
    abs.(signals) ./ matrices.volumes
end

# Across the axons, the water in the axons barely attenuates: the axons are too small for
# the water to move far in the gradient direction. Along the axons, the diffusion is free,
# and all compartments attenuate. The magnetization across the axons shows the restriction:

plot_field(mesh, magnetizations.across; decorations = false)

# The ADCs of the compartments follow from the homogenized ADC model, without solving the
# Bloch-Torrey equation for several b-values:

gradient = ScalarGradient(directions.across, profile; bvalue = 1000, gamma = model.gamma)
adcs = solve_hadc(model, matrices, gradient)
