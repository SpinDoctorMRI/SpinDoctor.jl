# # Solve BTPDE
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
using Random

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

# Set random seed for reproducibility
Random.seed!(123)

# Floating point type
T = Float64

# The built in geometry recipes allow for making various cell configurations.
# We here consider the case of three twisted axons immersed in an
# extracellular space (ECS). The axons are extruded from a 2D ground setup
# consisting of three random disks.

setup = CylinderSetup{T}(;
    ncell = 3,
    ecs_shape = :convex_hull,
    height = 40.0,
    bend = 0.0,
    twist = π / 4,
)

# We also define parameters for the cell compartments and ECS. The required
# parameters are:
#
# - 3D Diffusion tensors `D` (we here use an isotropic coefficient)
# - T2-relaxation times `T₂`
# - Initial spin densities `ρ`
# - Permeabilities `κ`
#
# The parameters `D`, `T₂`, and `ρ` are defined for each compartment, while
# `κ` is defined for each interface and outer boundary. The `DiskSetup` allows
# for multiple layers in each cell (here `nlayer = 1`), so `cell` and
# `cell_boundaries` are arrays of length `nlayer = 1`, while `cell_interfaces`
# is of length `nlayer - 1 = 0`.

coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3)], ecs = 0.002 * I(3)),
    T₂ = (; cell = [Inf], out = Inf, ecs = Inf),
    ρ = (; cell = [1.0], out = 1.0, ecs = 1.0),
    κ = (; cell_interfaces = zeros(0), cell_boundaries = [0.0], cell_ecs = 1e-4, ecs = 0.0),
    γ = 2.67513e-4,
)

# The following line creates a random cell configuration for our cylinder
# setup, generates a surface triangulation and calls TetGen to create a
# tetrahedral finite element mesh. The compartments and boundaries will be
# ordered in the same way as `coeffs`.

mesh, surfaces, cells = create_geometry(setup);

# The resulting mesh can be plotted in 3D provided a Makie backend is loaded.

plot_mesh(mesh)

# The mesh looks good, so we can proceed with the assembly our biological
# model and the associated finite element matrices.

model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)

# The Bloch-Torrey PDE takes a magnetic field gradient pulse sequence as an input. Here
# we consider a [`ScalarGradient`](@ref) with a [`PGSE`](@ref) time profile.

dir = [1.0, 0.0, 0.0]
profile = PGSE(2000.0, 6000.0)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
gradient = ScalarGradient(dir, profile, g)

# SpinDoctor provides a [`solve`](@ref) function, which has the same base signature for all
# diffusion MRI problems. The BTPDE is one such problem. They generally take a gradient
# sequence as an input.

btpde = BTPDE(; model, matrices)
ξ = solve(btpde, gradient)

# Here, `ξ` is a vector containing the complex-valued magnetization at all degrees of freedom
# at the echo time `TE`. We may compute the resulting signal as follows:

compute_signal(matrices.M, ξ)

# The global mass matrix `M` is used to compute the integral. We may however be interested in
# the compartment-wise signals. This requires splitting the magnetization field into the
# respective compartments. The compartment mass matrices are also available.

ξ_cmpts = split_field(mesh, ξ)
compute_signal.(matrices.M_cmpts, ξ_cmpts)

# The final magnetization can be visualized using the [`plot_field`](@ref) function.

plot_field(mesh, ξ)

# In this example, we have computed the complex transverse water proton magnetization field
# using the finite element method. The measured diffusion MRI signal is the integral of this
# field, and other quantities of interest, such as the apparent diffusion coefficient (ADC),
# or the effective diffusion tensor, may easily be obtained from this reference field.
# Directly solving the BTPDE is thus considered to be the "gold standard" for computing these
# quantities, as arbitrary precision may be obtained.

# However, this is also often the most computationally expensive approach. In the following
# examples, we will consider some other specialized methods provided by SpinDoctor, each
# having their own domains of validity, use cases, and computational footprints.
