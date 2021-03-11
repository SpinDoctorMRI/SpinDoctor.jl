"Diffusion MRI modelling"
module SpinDoctor

using Printf
using Parameters
using LinearAlgebra
using Arpack: eigs
using SparseArrays
using Statistics: mean

using DifferentialEquations
using Expokit: expmv!
using QuadGK
using Polyhedra: polyhedron, vrep, removevredundancy!, DefaultLibrary
using QHull: chull
import GLPK: Optimizer
using Triangle: constrained_triangulation
using TetGen
using WriteVTK
using Polynomials: fit

# using ForwardDiff
# can_dual(::Type{ComplexF64}) = true

# Setups
include("sequences.jl")
include("setup.jl")
include("prepare_pde.jl")

# Geometry
include("geometry/create_cells.jl")
include("geometry/create_surface_triangulation.jl")
include("geometry/create_mesh.jl")
include("geometry/deform_domain.jl")
include("geometry/split_mesh.jl")
include("geometry/create_fibonacci_sphere.jl")
include("geometry/create_directions.jl")
include("geometry/get_volumes.jl")

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


export integral, bvalue_no_q, echotime, PGSE, CosOGSE, SinOGSE, DoublePGSE, echotime,
    CellSetup, DomainSetup, ExperimentSetup, BTPDE, HADC, MF,
    prepare_pde,
    create_cells,
    create_surface_triangulation,
    create_mesh,
    deform_domain!,
    split_mesh,
    create_fibonacci_sphere,
    create_directions,
    get_cmpt_volumes,
    compute_laplace_eig, length2eig, eig2length,
    solve_btpde,
    solve_mf,
    savefield, savefield_time, save_btpde_results,
    fit_adc,
    Trapezoid, ImplicitEuler, ABDF2, QNDF, QNDF1, QNDF2, QBDF, QBDF1, QBDF2

end
