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
using SparseArrays
using Statistics: mean
using TetGen
using Triangle: constrained_triangulation
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
export Setup
export prepare_pde!
export prepare_experiments!
export create_geometry
export split_mesh
export create_directions
export get_cmpt_volumes
export compute_laplace_eig
export length2eig
export eig2length
export solve_btpde
export solve_mf
export savefield
export savefield_time
export save_btpde_results
export fit_adc
export Trapezoid
export ImplicitEuler
export ABDF2, QBDF, QBDF1, QBDF2, QNDF, QNDF1, QNDF2

# Setups
include("sequences.jl")
include("setup.jl")
include("prepare_pde.jl")

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
include("matrix_assembly/assemble_flux_matrix.jl")
include("matrix_assembly/assemble_flux_matrix_cmpt.jl")
include("matrix_assembly/couple_flux_matrix.jl")

# Eigendecomposition
include("eig2length.jl")
include("compute_laplace_eig.jl")

# Solvers
include("solve_btpde.jl")
include("solve_mf.jl")

# Postprocessing
include("fit_adc.jl")
include("savefield.jl")

end
