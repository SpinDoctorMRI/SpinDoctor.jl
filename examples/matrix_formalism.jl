# # Matrix Formalism
#
# In this example we will consider the matrix formalism approach for a geometry of cylinders.

# LSP indexing solution                                                          #src
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983 #src
if isdefined(@__MODULE__, :LanguageServer)                                       #src
    include("../src/SpinDoctor.jl")                                              #src
    using .SpinDoctor                                                            #src
end                                                                              #src

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
    nsidewall = 15,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    height = 30.0,
    bend = 0.0,
    twist = π / 6,
    ecs_shape = :convex_hull,
    ecs_ratio = 0.5,
)

# We also define coefficients for the different cell compartments `:in` (axon), `:out`
# (myelin), and `:ecs` (ECS).

coeffs = coefficients(
    setup;
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)

# We then proceed to build the geometry and finite element mesh.

mesh, = create_geometry(setup; recreate = true)
plot_mesh(mesh)

# The mesh looks good, so we may then proceed to assemble the biological model and the
# associated finite element matrices.

model = Model(; mesh, coeffs...);
matrices = assemble_matrices(model);

# We may also compute some useful quantities, including a scalar diffusion coefficient from
# the diffusion tensors.

volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)

# The eigenfunctions of the diffusive part of the Bloch-Torrey operator forms a good basis
# for the finite element function space. The basis may be truncated at a certain level, thus
# reducing the number of degrees of freedom. We here opt for 400 eigenfunctions.

laplace = Laplace(; model, matrices, neig_max = 400)
lap_eig = solve(laplace)

# The resulting eigenvalues may be represented as length scales, describing the wavelength
# of the eigenfunctions.

length_scales = eig2length.(lap_eig.values, D_avg)

# We may also further truncate the eigenfunction basis, if we are satisfied skipping
# features below a threshold length scale of 3 micrometers.

length_scale = 3
λ_max = length2eig(length_scale, D_avg)
lap_eig = limit_lengthscale(lap_eig, λ_max)

# Each of the resulting eigenfunctions is represented in the same way as the initial
# magnetization field `ρ`.

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

# We observe that the first functions have large features, while the higher-index functions
# have more rapidly varying features. We may now choose a gradient and compute the
# projection of magnetization field onto the truncated basis.

dir = [1.0, 0.0, 0.0]
profile = CosOGSE(5000.0, 5000.0, 2)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
gradient = ScalarGradient(dir, profile, g)

# The matrix formalism problem is solved in the same way as the [`BTPDE`](@ref). The time
# profile is approximated on 500 points, since it is non-constant.

mf = MatrixFormalism(; model, matrices, lap_eig)
ξ = solve(mf, gradient; ninterval = 500)

# The resulting magnetization field may be plotted.

plot_field(model.mesh, ξ)
