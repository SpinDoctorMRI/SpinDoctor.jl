using SpinDoctor
using LinearAlgebra

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

setup = SphereSetup(;
    name = "gaze-into-the-orb",
    ncell = 1,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    include_in = false,
    in_ratio = 0.6,
    ecs_shape = :no_ecs,
    ecs_ratio = 0.5,
)
coeffs = coefficients(
    setup;
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)

mesh, = create_geometry(setup; recreate = true)
plot_mesh(mesh)

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)

# Magnetic field gradient
TE = 5000.0
φ = -π / 6
R = [cos(φ) sin(φ) 0; -sin(φ) cos(φ) 0; 0 0 1]
g⃗(t) = 1.0 * R * [sin(2π * t / TE), sin(20π * t / TE) / 5, cos(2π * t / TE)]
gradient = GeneralGradient(; g⃗, TE)

plotter = Plotter{Float64}()

btpde = BTPDE(; model, matrices)
ξ = solve(btpde, gradient; callbacks = [plotter])
plotter.fig

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

