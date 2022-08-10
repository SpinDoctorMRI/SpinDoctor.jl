# # Custom gradients
#
# Here we will consider a custom gradient applied to a spherical geometry.
#
# We start by loading SpinDoctor and a Makie plotting backend.

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

T = Float64
setup = SphereSetup{T}(; ncell = 1, layersizes = [0.6, 1.0], rmin = 5.0, rmax = 5.0)
nlayer = length(setup.layersizes)
coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3) for _ = 1:nlayer], ecs = 0.002 * I(3)),
    T₂ = (; cell = fill(Inf, nlayer), ecs = Inf),
    ρ = (; cell = fill(1.0, nlayer), ecs = 1.0),
    κ = (;
        cell_interfaces = fill(1e-4, nlayer - 1),
        cell_boundaries = fill(0, nlayer),
        cell_ecs = 1e-4,
        ecs = 0,
    ),
    γ = 2.67513e-4,
)

# We then proceed to build the geometry and finite element mesh.

mesh, = create_geometry(setup)
plot_mesh(mesh)

# The mesh looks good, so we may then proceed to assemble the biological model and the
# associated finite element matrices.

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)

# The Bloch-Torrey PDE takes a magnetic field gradient pulse sequence as an input. We may
# define our custom three dimensional gradient sequence (given in T/m) as a simple Julia
# function. The echo time `TE` (given in microseconds) is also needed.

## Magnetic field gradient
TE = 5000.0
φ = -π / 6
R = [cos(φ) sin(φ) 0; -sin(φ) cos(φ) 0; 0 0 1]
gvec(t) = 1.0 * R * [sin(2π * t / TE), sin(20π * t / TE) / 5, cos(2π * t / TE)]
gradient = GeneralGradient(; gvec, TE)

# In order to follow the evolution of the solution during time stepping, we add a
# [`Plotter`](@ref) to a list of callbacks. Other available callbacks are [`Printer`](@ref)
# for showing time stepping information, and [`VTKWriter`](@ref) for saving the solution time
# series for visualization in ParaView.

plotter = Plotter{Float64,3}()

# We may then define the problem and solve for our gradient (with the callback).

btpde = BTPDE(; model, matrices)
ξ = solve(btpde, gradient; callbacks = [plotter])
plotter.fig

# [`MatrixFormalism`](@ref) also accepts [`GeneralGradient`](@ref)s. However, the ADC methods
# currently only support [`ScalarGradient`](@ref)s.
