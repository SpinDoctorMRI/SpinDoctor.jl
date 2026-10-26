# # [Validation against analytical solutions](@id validation)
#
# For layered plates, cylinders and spheres, the matrix formalism has analytical
# eigenfunctions [Grebenkov2010a](@cite), which give reference signals for the finite element
# solutions. Here, we compare the signals of four permeable plates.

using SpinDoctor
using LinearAlgebra
using CairoMakie

# Four plates of 5 µm in 2D, with permeable interfaces. With a gradient across the plates,
# the problem is one-dimensional.

setup = PlateSetup(; depth = 20.0, widths = fill(5.0, 4), refinement = 0.5)
coeffs = coefficients(
    setup;
    D = fill(0.002 * I(2), 4),
    T2 = fill(Inf, 4),
    rho = fill(1.0, 4),
    kappa = (; interfaces = fill(1e-4, 3), boundaries = fill(0.0, 4)),
    gamma = 2.67513e-4,
)
mesh, = create_geometry(setup)
plot_mesh(mesh)

# The analytical eigenvalues and eigenfunctions of the layered plates, with length scales
# above 0.3 µm:

volumes = compartment_volumes(mesh)
laplace = solve_analytical_laplace(;
    analytical_coefficients(setup, coeffs)...,
    eiglim = length2eig(0.3, 0.002),
    eigstep = 1e-8,
)
analytical(gradient) = real(solve_analytical_mf(laplace, volumes, gradient)) / sum(volumes);

# The finite element signals, for meshes with smaller and smaller triangles:

profile = PGSE(5000.0, 15000.0)
bvalues = 0:500:5000
gradients =
    [ScalarGradient([1.0, 0.0], profile; bvalue = b, gamma = coeffs.gamma) for b in bvalues]
refinements = [2.0, 0.5, 0.125]
fem = map(refinements) do refinement
    setup = PlateSetup(; depth = 20.0, widths = fill(5.0, 4), refinement)
    mesh, = create_geometry(setup)
    model = Model(; mesh, coeffs...)
    matrices = assemble_matrices(model)
    solver = ThetaSolver(; timestep = 10.0)
    signals = map(gradients) do gradient
        u = solve_btpde(model, matrices, gradient, solver)
        real(compute_signal(matrices.M, u)) / sum(matrices.volumes)
    end
    (; npoint = sum(p -> size(p, 2), mesh.points), signals)
end
reference = analytical.(gradients)

# The signals agree, and the error decreases with the mesh size, about as the square of the
# size of the triangles (second order):

fig = Figure(; size = (900, 380))
ax = Axis(fig[1, 1]; xlabel = "b [s/mm²]", ylabel = "S / S₀", yscale = log10)
lines!(ax, bvalues, reference; color = :black, label = "Analytical")
scatter!(ax, bvalues, fem[end].signals; label = "Finite elements")
axislegend(ax; position = :lb)
ax = Axis(fig[1, 2]; xlabel = "b [s/mm²]", ylabel = "Relative error", yscale = log10)
for (refinement, f) in zip(refinements, fem)
    error = abs.(f.signals .- reference)[2:end] ./ reference[2:end]
    scatterlines!(ax, bvalues[2:end], error; label = "$(f.npoint) points")
end
axislegend(ax; position = :rb)
fig
