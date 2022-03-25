# # High angular resolution diffusion imaging
#
# In this example we will compute the signal using the same gradient sequence in many
# different directions.

# ## Building a biological model
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

# Here we create a recipe for five stacked plates with isotropic diffusion tensors. They
# should allow for free diffusion in the horizontal direction, but a rather restricted
# vertical diffusion with the permeable membranes.

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

# We then proceed to build the geometry and finite element mesh.

mesh, = create_geometry(setup; recreate = true)
plot_mesh(mesh)

# The mesh looks good, so we may then proceed to assemble the biological model and the
# associated finite element matrices.

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model);

# We may also compute some useful quantities, including a scalar diffusion coefficient from
# the diffusion tensors.

volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)
ncompartment = length(model.mesh.points)

# The gradient pulse sequence will be a PGSE with both vertical and horizontal components.
# This allows for both restricted vertical diffusion and almost unrestricted horizontal
# diffusion. The different approaches should hopefully confirm this behaviour.

directions = unitsphere(200)
profile = PGSE(2500.0, 4000.0)
b = 1000
g = √(b / int_F²(profile)) / model.γ
gradients = [ScalarGradient(d, profile, g) for d ∈ eachcol(directions)]

# We may solve the BTPDE for each gradient.

btpde = BTPDE(; model, matrices)
solver = IntervalConstantSolver(; timestep = 10.0)
ξ, = solve_multigrad(btpde, gradients, solver)

# The signals are computed from the magnetization field through quadrature.

signal = [abs(compute_signal(matrices.M, ξ)) for ξ ∈ ξ] 

# We may plot the directionalized signal.

plot_hardi(directions, signal)

# The signal attenuates the most in the horizontal direction, as that is where diffusion is
# restricted the least.
