"Diffusion MRI modelling"
module SpinDoctor

using Arpack: eigs
using Expokit: expmv!
using GLMakie
using GLPK: Optimizer
using LinearAlgebra
using MiniQhull: delaunay
using OrdinaryDiffEq
using Parameters
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
using WriteVTK

export TimeProfile,
    PGSE,
    CosOGSE,
    SinOGSE,
    DoublePGSE,
    integral,
    bvalue_no_q,
    intervals,
    get_interval,
    constant_intervals,
    is_constant,
    echotime
export Setup, CylinderSetup, SphereSetup, NeuronSetup, Experiment
export get_coefficients
export create_cells
export create_surfaces_cylinder
export create_surfaces_sphere
export create_surfaces_neuron
export create_geometry
export create_model
export split_mesh
export create_directions
export get_cmpt_volumes, get_mesh_volumes
export compute_laplace_eig
export length2eig
export eig2length
export Model
export savefield
export savefield_time
export save_btpde_results
export solve_analytical
export solve_btpde
export solve_btpde_midpoint
export solve_hadc
export solve_mf
export solve_karger
export fit_adc
export fit_tensors
export plot_mesh
export get_values

# Rexport default ODE solver
export QNDF, MagnusGL6

# Datatypes
include("datatypes/sequences.jl")
include("datatypes/experiment.jl")
include("datatypes/femesh.jl")
include("datatypes/gradient.jl")
include("datatypes/model.jl")
include("datatypes/create_model.jl")
include("datatypes/setup.jl")
include("datatypes/get_coefficients.jl")

# Geometry
include("geometry/call_tetgen.jl")
include("geometry/convexhull.jl")
include("geometry/create_cells.jl")
include("geometry/create_directions.jl")
include("geometry/create_fibonacci_sphere.jl")
include("geometry/create_geometry.jl")
include("geometry/create_surfaces_cylinder.jl")
include("geometry/create_surfaces_neuron.jl")
include("geometry/create_surfaces_sphere.jl")
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

# Matrix assembly
include("matrix_assembly/assemble_mass_matrix.jl")
include("matrix_assembly/assemble_stiffness_matrix.jl")
include("matrix_assembly/assemble_flux_matrices.jl")
include("matrix_assembly/assemble_flux_matrix.jl")
include("matrix_assembly/couple_flux_matrix.jl")

# Matrix formalism
include("matrix_formalism/eig2length.jl")
include("matrix_formalism/compute_laplace_eig.jl")
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
include("analytical/solve_analytical.jl")

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

# Utils
include("utils/get_values.jl")

end
