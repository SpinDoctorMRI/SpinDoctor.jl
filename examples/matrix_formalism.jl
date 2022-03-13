# # Matrix Formalism
#
# In this example we will consider the matrix formalism approach for a geometry of cylinders.

using SpinDoctor
using LinearAlgebra

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

# LSP indexing solution                                                          #src
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983 #src
if isdefined(@__MODULE__, :LanguageServer)                                       #src
    include("../src/SpinDoctor.jl")                                              #src
    using .SpinDoctor                                                            #src
end                                                                              #src

setup = CylinderSetup(;
    name = "Slice",
    ncell = 10,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    height = 1.0,
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

mesh, = create_geometry(setup; recreate = true);
model = Model(; mesh, coeffs...);
matrices = assemble_matrices(model);

volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)

plot_mesh(model.mesh)

laplace = Laplace{T}(; model, matrices, neig_max = 400)
lap_eig = solve(laplace)
length_scales = eig2length.(lap_eig.values, D_avg)

length_scale = 3
λ_max = length2eig(length_scale, D_avg)
lap_eig = limit_lengthscale(lap_eig, λ_max)

fig = Figure()
for i = 1:3, j = 1:6
    icmpt = 6(i - 1) + j
    ϕ = lap_eig.funcs[icmpt]
    ax = Axis(fig[i, j])
    # TODO: Plot different lap_eig
end
fig

dir = [1.0, 0.0, 0.0]
profile = CosOGSE(5000.0, 5000.0, 2)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
gradient = ScalarGradient(dir, profile, g)

mf = MatrixFormalism(; model, matrices, lap_eig)
ξ = solve(mf, gradient; ninterval = 500)

plot_field(model.mesh, ξ)
