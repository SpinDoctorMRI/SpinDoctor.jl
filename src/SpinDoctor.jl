"""
    SpinDoctor

Diffusion MRI simulation toolbox.
"""
module SpinDoctor

using DiffEqCallbacks: FunctionCallingCallback
using Expokit: expmv!
using KrylovKit: KrylovKit, InnerProductVec
using LinearAlgebra
using OrdinaryDiffEqBDF: QNDF
using Parameters
using Polynomials: fit
using Printf
using QuadGK
using Random
using Roots: find_zeros
using SciMLBase: SciMLBase, AbstractODEAlgorithm, ODEFunction, ODEProblem
using SparseArrays
using SpecialFunctions
using StaticArrays
using Statistics: mean
using TetGen: RawTetGenIO, facetlist!, tetrahedralize
using Triangulate: TriangulateIO, triangulate
using WriteVTK

# Gradients
export PGSE, CosOGSE, SinOGSE, DoublePGSE, GeneralGradient, ScalarGradient
export integral, int_F², intervals, isconstant, echotime

export assemble_matrices, assemble_mass_matrix, assemble_stiffness_matrix
export get_mesh_volumes, get_mesh_surface, get_mesh_surfacenormals
export split_mesh
export split_field
export get_cmpt_volumes
export compartment_mesh, box_faces

export Model
export initial_conditions

export savefield

export compute_adc_sta
export fit_adc
export fit_tensors
export compartment_weights
export plot_surfaces, plot_mesh, plot_field, plot_hardi

export solve_btpde, solve_hadc, solve_karger, solve_laplace, solve_mf
export solve_analytical_laplace, solve_analytical_mf
export QNDF, IntervalConstantSolver

# Matrix formalism
export length2eig, eig2length, limit_lengthscale, compute_mf_diffusion_tensor

# Utils
export unitcircle, unitsphere
export compute_signal

# Callbacks
export Printer, VTKWriter, Plotter

# Recipes
export NoECS, BoxECS, ConvexHullECS
export AbstractSetup,
    PlateSetup,
    ExtrusionSetup,
    DiskSetup,
    SphereSetup,
    NeuronSetup,
    CylinderSetup,
    SlabSetup
export coefficients, analytical_coefficients
export create_geometry, save_geometry, load_geometry

include("directions.jl")
include("gradients.jl")
include("gradient_integrals.jl")

include("mesh.jl")
include("meshing.jl")
include("io.jl")
include("model.jl")

include("assembly.jl")
include("matrices.jl")

include("btpde.jl")
include("karger.jl")
include("adc.jl")
include("matrix_formalism.jl")
include("analytical.jl")
include("postprocess.jl")
include("callbacks.jl")

include("recipes/setups.jl")
include("recipes/coefficients.jl")
include("recipes/cells.jl")
include("recipes/surfaces.jl")
include("recipes/geometry.jl")

include("extensions.jl")

end
