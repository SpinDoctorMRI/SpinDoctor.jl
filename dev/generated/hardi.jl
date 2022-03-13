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

directions = unitsphere(200)
profile = PGSE(2500.0, 4000.0)
b = 1000
g = √(b / int_F²(profile)) / model.γ
gradients = [ScalarGradient(d, profile, g) for d ∈ eachcol(directions)]

btpde = BTPDE(; model, matrices)
solver = IntervalConstantSolver(; timestep = 10.0)
ξ, = solve_multigrad(btpde, gradients, solver)

signal = [abs(compute_signal(matrices.M, ξ)) for ξ ∈ ξ]

plot_hardi(directions, signal)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

