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
using Printf

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

setup = DiskSetup(;
    ncell = 5,
    nsidewall = 20,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    layersizes = [0.3, 0.6, 1.0],
    ecs = ConvexHullECS(; margin = 2.0),
    refinement = 0.2,
)

# We also define coefficients for the three cell compartments and the ECS.

nlayer = length(setup.layersizes)
coeffs = coefficients(
    setup;
    D = (; cell = [0.001 * I(2) for _ = 1:nlayer], ecs = 0.002 * I(2)),
    T₂ = (; cell = fill(Inf, nlayer), ecs = Inf),
    ρ = (; cell = fill(1.0, nlayer), ecs = 1.0),
    κ = (;
        cell_interfaces = fill(1e-4, nlayer - 1),
        cell_boundaries = fill(0.0, nlayer),
        cell_ecs = 1e-5,
        ecs = 0.0,
    ),
    γ = 2.67513e-4,
)

# We then proceed to build the geometry and finite element mesh.

mesh, surfaces, cells = create_geometry(setup; recreate = true)
plot_mesh(mesh, 1:3)
plot_mesh(mesh)

# The mesh looks good, so we may then proceed to assemble the biological model
# and the associated finite element matrices.

model = Model(; mesh, coeffs...);
matrices = assemble_matrices(model);

# We may also compute some useful quantities, including a scalar diffusion
# coefficient from the diffusion tensors.

volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 2 * tr.(model.D)' * volumes / sum(volumes)

# The eigenfunctions of the diffusive part of the Bloch-Torrey operator forms
# a good basis for the finite element function space. The basis may be
# truncated at a certain level, thus reducing the number of degrees of
# freedom. We here opt for 400 eigenfunctions.

laplace = Laplace(; model, matrices, neig_max = 200)
lap_eig = solve(laplace)

# The resulting eigenvalues may be represented as length scales, describing
# the wavelength of the eigenfunctions.

length_scales = eig2length.(lap_eig.values, D_avg)

# We may also further truncate the eigenfunction basis, if we are satisfied
# skipping features below a threshold length scale of 3 micrometers.

length_scale = 3
λ_max = length2eig(length_scale, D_avg)
lap_eig = limit_lengthscale(lap_eig, λ_max)

# Each of the resulting eigenfunctions is represented in the same way as the initial
# magnetization field `ρ`.

ncompartment = length(mesh.points)
fig = Figure()
for i = 1:3, j = 1:4
    ieig = 1 + 6(i - 1) + j
    ϕ_cmpts = split_field(mesh, lap_eig.funcs[:, ieig])
    ax = Axis(fig[i, j]; title = @sprintf("n = %d, ℓ = %.1f", ieig, length_scales[ieig]))
    for icmpt ∈ 1:ncompartment
        f = mesh.elements[icmpt]
        p = mesh.points[icmpt]
        mesh!(ax, p', f'; color = ϕ_cmpts[icmpt], shading = false)
    end
end
fig

# We observe that the first functions have large features, while the higher-index functions
# have more rapidly varying features. We may now choose a gradient and compute the
# projection of magnetization field onto the truncated basis.

dir = [1.0, 0.0]
profile = CosOGSE(5000.0, 5000.0, 2)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
gradient = ScalarGradient(dir, profile, g)

# The matrix formalism problem is solved in the same way as the
# [`BTPDE`](@ref). The time profile is approximated on 500 points, since it is
# non-constant.

mf = MatrixFormalism(; model, matrices, lap_eig)
ξ = solve(mf, gradient; ninterval = 500)

# The resulting magnetization field may be plotted.

plot_field(model.mesh, ξ)
