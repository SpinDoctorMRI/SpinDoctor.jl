# # [Matrix formalism](@id mf-example)
#
# The matrix formalism [Agdestein2021](@cite) expands the magnetization in the Laplace
# eigenfunctions of the geometry (with the permeable membranes). The eigenfunctions with
# long length scales describe the signal, and truncating the basis at a length scale gives
# a small model, which is solved quickly for many gradient sequences.

using SpinDoctor
using LinearAlgebra
using Random
using Printf
using CairoMakie

# ## Model
#
# Five cells with three layers (radii of 30 %, 60 % and 100 % of the cell radius) in an
# ECS, in 2D:

ncell = 5
setup = DiskSetup(;
    ncell,
    rmin = 2.0,
    rmax = 6.0,
    nsidewall = 20,
    dmin = 0.2,
    dmax = 0.3,
    layersizes = [0.3, 0.6, 1.0],
    ecs = ConvexHullECS(; margin = 2.0),
    refinement = 0.2,
)
mesh, = create_geometry(setup; rng = Xoshiro(2))
layercolors = [:crimson, :orange, :steelblue]
plot_mesh(mesh; colors = [repeat(layercolors; inner = ncell); :lightgray])

#-

nlayer = 3
coeffs = coefficients(
    setup;
    D = (; cell = [0.001 * I(2) for _ = 1:nlayer], ecs = 0.002 * I(2)),
    T2 = (; cell = fill(Inf, nlayer), ecs = Inf),
    rho = (; cell = fill(1.0, nlayer), ecs = 1.0),
    kappa = (;
        cell_interfaces = fill(1e-4, nlayer - 1),
        cell_boundaries = fill(0.0, nlayer),
        cell_ecs = 1e-5,
        ecs = 0.0,
    ),
    gamma = 2.67513e-4,
)
model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)
D_mean = sum(matrices.volumes .* tr.(model.D) ./ 2) / sum(matrices.volumes)

# ## Laplace eigenfunctions
#
# The 400 smallest eigenvalues, and their length scales (the wavelengths of the
# eigenfunctions):

laplace = solve_laplace(model, matrices; neig_max = 400)
lengths = eig2length.(laplace.values, D_mean)

# The first eigenfunctions vary slowly, and the next ones more and more rapidly:

indices = [2, 3, 5, 10, 15, 20, 30, 40, 50, 100, 200, 400]
fig = Figure(; size = (1000, 650))
for (k, n) in enumerate(indices)
    ax = Axis(
        fig[div(k-1, 4)+1, mod1(k, 4)];
        title = @sprintf("n = %d, ℓ = %.1f µm", n, lengths[n]),
        aspect = DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)
    plot_field!(ax, mesh, laplace.funcs[:, n]; part = :real)
end
fig

# ## Signal
#
# The reference signal of the Bloch-Torrey equation for an oscillating gradient:

profile = CosOGSE(5000.0, 5000.0, 2)
gradient = ScalarGradient([1.0, 0.0], profile; bvalue = 1000, gamma = model.gamma)
u = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 10.0))
signal = compute_signal(matrices.M, u)

# The matrix formalism signal with the eigenfunctions of length scales above `ℓ`, for
# decreasing `ℓ`:

lmins = [20.0, 10.0, 6.0, 4.0, 3.0, 2.0, 1.5]
mf = map(lmins) do lmin
    truncated = limit_lengthscale(laplace, length2eig(lmin, D_mean))
    u_mf = solve_mf(model, matrices, truncated, gradient)
    (; neig = length(truncated.values), signal = compute_signal(matrices.M, u_mf))
end
errors = [abs(m.signal - signal) / abs(signal) for m in mf]

# The error decreases with the number of eigenfunctions, that is with the resolved length
# scale:

fig = Figure(; size = (600, 380))
ax = Axis(
    fig[1, 1];
    xlabel = "Number of eigenfunctions",
    ylabel = "Relative error of the signal",
    xscale = log10,
    yscale = log10,
    xautolimitmargin = (0.05, 0.2),
)
scatterlines!(ax, [m.neig for m in mf], errors)
for (m, lmin, error) in zip(mf, lmins, errors)
    text!(ax, m.neig, error; text = " ℓ = $lmin µm", fontsize = 12)
end
fig

# The magnetization of the matrix formalism, with the eigenfunctions above 2 µm:

truncated = limit_lengthscale(laplace, length2eig(2.0, D_mean))
plot_field(mesh, solve_mf(model, matrices, truncated, gradient))
