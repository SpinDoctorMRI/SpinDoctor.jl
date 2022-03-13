using SpinDoctor
using LinearAlgebra

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

setup = CylinderSetup(;
    name = "some-very-real-axons",
    ncell = 3,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    height = 40.0,
    bend = 0.0,
    twist = π / 4,
    include_in = false,
    in_ratio = 0.6,
    ecs_shape = :convex_hull,
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

mesh, = create_geometry(setup)

plot_mesh(mesh)

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)

dir = [1.0, 0.0, 0.0]
profile = PGSE(2000.0, 6000.0)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
gradient = ScalarGradient(dir, profile, g)

btpde = BTPDE(; model, matrices)
ξ = solve(btpde, gradient)

compute_signal(matrices.M, ξ)

ξ_cmpts = split_field(mesh, ξ)
compute_signal.(matrices.M_cmpts, ξ_cmpts)

plot_field(mesh, ξ)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

