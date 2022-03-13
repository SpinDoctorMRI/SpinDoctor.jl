using SpinDoctor
using LinearAlgebra

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

setup = CylinderSetup(;
    name = "Slice",
    ncell = 3,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    height = 1.0,
    bend = 0.0,
    twist = 0.0,
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

mesh, = create_geometry(setup; recreate = true)
plot_mesh(mesh)

model = Model(; mesh, coeffs...);
matrices = assemble_matrices(model);

volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)

laplace = Laplace(; model, matrices, neig_max = 400)
lap_eig = solve(laplace)

length_scales = eig2length.(lap_eig.values, D_avg)

length_scale = 3
λ_max = length2eig(length_scale, D_avg)
lap_eig = limit_lengthscale(lap_eig, λ_max)

ncompartment, nboundary = size(mesh.facets)
fig = Figure()
for i = 1:3, j = 1:4
    ieig = 6(i - 1) + j
    ϕ_cmpts = split_field(mesh, lap_eig.funcs[:, ieig])
    ax = Axis3(fig[i, j]; title = "n = $ieig, ℓ = $(length_scales[ieig])", aspect = :data)
    nboundary = size(mesh.facets, 2)
    scene = nothing
    first = true
    for icmpt = 1:ncompartment, iboundary = 1:nboundary
        facets = mesh.facets[icmpt, iboundary]
        points = mesh.points[icmpt]
        mesh!(ax, points', facets', color = ϕ_cmpts[icmpt], shading = false)
    end
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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

