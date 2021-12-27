"Diffusion MRI modelling"
module SpinDoctor

using Arpack: eigs
using DiffEqCallbacks: PresetTimeCallback, FunctionCallingCallback
using Expokit: expmv!
using GLPK: Optimizer
using LinearAlgebra
using Makie
using MiniQhull: delaunay
using OrdinaryDiffEq: OrdinaryDiffEq
using OrdinaryDiffEq: DiffEqArrayOperator, DiscreteCallback, ODEProblem, ODEFunction, QNDF, MagnusGL6
using Polynomials: fit
using Printf
using QuadGK
using Roots: find_zeros
using SparseArrays
using SpecialFunctions
using StaticArrays
using Statistics: mean
using TetGen: RawTetGenIO, facetlist!, tetrahedralize
using Triangulate: TriangulateIO, triangulate
using UnPack
using WriteVTK

# Gradidents
export PGSE, CosOGSE, SinOGSE, DoublePGSE, GeneralGradient, ScalarGradient
export integral, int_F², intervals, isconstant, echotime

export assemble_matrices
export create_cells
export create_model
export split_mesh
export split_field
export get_cmpt_volumes, get_mesh_volumes
export length2eig, eig2length
export Model
export savefield
export solve, solve_multigrad
export fit_adc
export fit_tensors
export plot_mesh, plot_field
export compute_lap_eig_analytical
export limit_lengthscale

export IntervalConstanBTPDE, GeneralBTPDE, HADC, Karger
export Laplace, MatrixFormalism, AnalyticalLaplace, AnalyticalMatrixFormalism
export solve

# Utils
export unitcircle, unitsphere
export compute_signal

# Callbacks
export VTKWriter, Plotter
export initialize!, update!, finalize!

# Recipes
export AbstractSetup, CylinderSetup, SphereSetup, NeuronSetup
export coefficients, analytical_coefficients, radial_dimension
export create_surfaces_cylinder, create_surfaces_neuron, create_surfaces_sphere, create_geometry

# Rexport default ODE solver
export QNDF, MagnusGL6

# Magnetic field gradient sequences
include("gradients/sequences.jl")
include("gradients/gradient.jl")
include("gradients/integral.jl")
include("gradients/int_F2.jl")
include("gradients/intervals.jl")
include("gradients/echotime.jl")
include("gradients/isconstant.jl")

# Datatypes
include("datatypes/femesh.jl")
include("datatypes/model.jl")

# Problems
include("problems/problems.jl")
include("problems/output_type.jl")
include("problems/solve.jl")
include("problems/solve_multigrad.jl")

# Geometry
include("geometry/call_tetgen.jl")
include("geometry/convexhull.jl")
include("geometry/create_cells.jl")
include("geometry/create_fibonacci_sphere.jl")
include("geometry/deform_domain.jl")
include("geometry/get_volumes.jl")
include("geometry/get_mesh_surface.jl")
include("geometry/get_mesh_surfacenormals.jl")
include("geometry/gmesh2fem.jl")
include("geometry/read_cells.jl")
include("geometry/read_surfaces.jl")
include("geometry/read_tetgen.jl")
include("geometry/save_cells.jl")
include("geometry/save_surfaces.jl")
include("geometry/save_tetgen.jl")
include("geometry/split_mesh.jl")
include("geometry/split_field.jl")

# Matrix assembly
include("matrix_assembly/assemble_mass_matrix.jl")
include("matrix_assembly/assemble_stiffness_matrix.jl")
include("matrix_assembly/assemble_flux_matrices.jl")
include("matrix_assembly/assemble_flux_matrix.jl")
include("matrix_assembly/couple_flux_matrix.jl")
include("matrix_assembly/assemble_matrices.jl")

# Matrix formalism
include("matrix_formalism/eig2length.jl")
include("matrix_formalism/limit_lengthscale.jl")
include("matrix_formalism/solve_laplace.jl")
include("matrix_formalism/solve_mf.jl")

# Analytical
include("analytical/alpha_func.jl")
include("analytical/compute_bc.jl")
include("analytical/compute_beta.jl")
include("analytical/compute_int_I.jl")
include("analytical/compute_int_J.jl")
include("analytical/compute_int_K.jl")
include("analytical/compute_JY.jl")
include("analytical/compute_nu.jl")
include("analytical/compute_v.jl")
include("analytical/find_alpha.jl")
include("analytical/solve_analytical_laplace.jl")
include("analytical/solve_analytical_matrix_formaliism.jl")

# BTPDE
include("btpde/solve_btpde.jl")
include("btpde/solve_btpde_midpoint.jl")
include("btpde/solve_karger.jl")

# ADC
include("adc/solve_hadc.jl")

# Postprocessing
include("postprocess/fit_adc.jl")
include("postprocess/fit_tensors.jl")
include("postprocess/savefield.jl")

# Plot
include("plot/plot_mesh.jl")
include("plot/plot_field.jl")

# Utils
include("utils/unitcircle.jl")
include("utils/unitsphere.jl")
include("utils/compute_signal.jl")

# Callbacks
include("callbacks/callbacks.jl")
include("callbacks/initialize.jl")
include("callbacks/update.jl")
include("callbacks/finalize.jl")

# Recipes
include("recipes/setup.jl")
include("recipes/coefficients.jl")
include("recipes/create_surfaces_cylinder.jl")
include("recipes/create_surfaces_neuron.jl")
include("recipes/create_surfaces_sphere.jl")
include("recipes/create_geometry.jl")
include("recipes/radial_dimension.jl")
include("recipes/analytical_coefficients.jl")

end
