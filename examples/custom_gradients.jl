# # [Custom gradient sequences](@id custom-gradients)
#
# Besides the scalar sequences (PGSE, double PGSE, OGSE), the gradient may be any function
# of time, with a varying direction ([`GeneralGradient`](@ref)). Here, we follow the
# magnetization of a cell with a nucleus during such a sequence.

using SpinDoctor
using LinearAlgebra
using CairoMakie

# A spherical cell with a nucleus (the inner 60 % of the radius), in an ECS:

setup = SphereSetup(; ncell = 1, layersizes = [0.6, 1.0], rmin = 5.0, rmax = 5.0)
nlayer = 2
coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3) for _ = 1:nlayer], ecs = 0.002 * I(3)),
    T2 = (; cell = fill(Inf, nlayer), ecs = Inf),
    rho = (; cell = fill(1.0, nlayer), ecs = 1.0),
    kappa = (;
        cell_interfaces = fill(1e-4, nlayer - 1),
        cell_boundaries = fill(0, nlayer),
        cell_ecs = 1e-4,
        ecs = 0,
    ),
    gamma = 2.67513e-4,
)
mesh, = create_geometry(setup)
plot_mesh(mesh; colors = [:crimson, (:steelblue, 0.5)])

#-

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model);

# The gradient (T/m) is a Julia function of time (µs), here a rotating field, with the echo
# time `TE`:

TE = 5000.0
phi = -pi / 12
R = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]
gvec(t) = 1.0 * R * [sin(2pi * t / TE), sin(20pi * t / TE) / 5, cos(2pi * t / TE)]
gradient = GeneralGradient(; gvec, TE)

# The [`Plotter`](@ref) callback follows the time stepping: the gradient, the signal
# attenuation, and the magnitude and phase of the magnetization. Other callbacks print
# information ([`Printer`](@ref)) or save the magnetization for ParaView
# ([`VTKWriter`](@ref)).

plotter = Plotter()
u = solve_btpde(model, matrices, gradient; callbacks = [plotter])
plotter.fig

# With an interactive backend (`using GLMakie`), the figure is shown and updated during
# the time stepping.
