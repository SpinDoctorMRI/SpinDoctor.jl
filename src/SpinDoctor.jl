"Diffusion MRI modelling"
module SpinDoctor

using Arpack: eigs
using DifferentialEquations
using Expokit: expmv!
using GLPK: Optimizer
using LinearAlgebra
using Parameters
using Polyhedra: DefaultLibrary, polyhedron, removevredundancy!, vrep # GenericLinearAlgebra problem
using Polynomials: fit
using Printf
using QHull: chull # GenericLinearAlgebra problem
using QuadGK
using Roots: find_zeros
using SparseArrays
using SpecialFunctions
using StaticArrays
using Statistics: mean
using TetGen
using Triangle: constrained_triangulation
using UnicodePlots
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
export solve_mf
export fit_adc
export Trapezoid
export ImplicitEuler
export QNDF

# Setups
include("datatypes/sequences.jl")
include("datatypes/experiment.jl")
include("datatypes/femesh.jl")
include("datatypes/gradient.jl")
include("datatypes/model.jl")
include("datatypes/create_model.jl")
include("datatypes/setup.jl")
include("get_coefficients.jl")

# Geometry
include("geometry/call_tetgen.jl")
include("geometry/create_cells.jl")
include("geometry/create_directions.jl")
include("geometry/create_fibonacci_sphere.jl")
include("geometry/create_geometry.jl")
include("geometry/create_surfaces_cylinder.jl")
include("geometry/create_surfaces_neuron.jl")
include("geometry/create_surfaces_sphere.jl")
include("geometry/deform_domain.jl")
include("geometry/get_volumes.jl")
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

# Eigendecomposition
include("eig2length.jl")
include("compute_laplace_eig.jl")

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

# Solvers
include("solve_btpde.jl")
include("solve_btpde_midpoint.jl")
include("solve_mf.jl")
include("solve_analytical.jl")

# Postprocessing
include("fit_adc.jl")
include("savefield.jl")

end
