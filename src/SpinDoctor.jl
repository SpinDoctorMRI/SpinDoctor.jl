"""
    SpinDoctor

Diffusion MRI simulation toolbox: finite element simulations of the Bloch-Torrey equation
and related models in geometries of cells and extracellular space.

The main functions and types are exported. Lower level building blocks (meshes, assembly,
geometry files, gradient sequence internals, callback hooks) are public but not exported,
and are used as `SpinDoctor.name`.
"""
module SpinDoctor

using DiffEqCallbacks: FunctionCallingCallback
using Expokit: expmv!
using KrylovKit: KrylovKit, InnerProductVec
using LinearAlgebra
using OrdinaryDiffEqBDF: QNDF
using Parameters: @with_kw
using Polynomials: fit
using Printf: @sprintf
using QuadGK: quadgk
using Random: Random
using Roots: find_zeros
using SciMLBase: SciMLBase, AbstractODEAlgorithm, ODEFunction, ODEProblem
using SparseArrays
using SpecialFunctions: besselj, bessely, erf, sphericalbesselj, sphericalbessely
using StaticArrays: SMatrix, SVector
using Statistics: mean
using TetGen: RawTetGenIO, facetlist!, tetrahedralize
using Triangulate: TriangulateIO, triangulate
using WriteVTK:
    MeshCell,
    VTKCellTypes,
    VTKFieldData,
    VTKPointData,
    paraview_collection,
    vtk_grid,
    vtk_multiblock,
    vtk_save

# Geometry
export PlateSetup, DiskSetup, SphereSetup, CylinderSetup, SlabSetup, ExtrusionSetup
export NeuronSetup, NoECS, BoxECS, ConvexHullECS
export create_geometry, save_geometry, load_geometry
export compartment_mesh, box_faces, compartment_volumes

# Model and matrices
export coefficients, analytical_coefficients, Model, initial_conditions
export assemble_matrices

# Gradients
export PGSE, DoublePGSE, CosOGSE, SinOGSE, ScalarGradient, GeneralGradient, bvalue
export unitcircle, unitsphere

# Solvers
export solve_btpde, IntervalConstantSolver, QNDF, solve_hadc, compute_adc_sta
export solve_karger, solve_laplace, solve_mf, solve_analytical_laplace, solve_analytical_mf
export eig2length, length2eig, limit_lengthscale, compute_mf_diffusion_tensor

# Post-processing and output
export compute_signal, split_field, compartment_weights, fit_adc, fit_tensors, save_field
export Printer, VTKWriter, Plotter
export plot_surfaces, plot_mesh, plot_field, plot_hardi

# Public, but not exported
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("""
        public FEMesh,
            split_mesh,
            facet_elements,
            element_geometry,
            facet_geometry,
            assemble_mass_matrix,
            assemble_stiffness_matrix,
            assemble_flux_matrix,
            AbstractSetup,
            AbstractECS,
            create_cells,
            create_surfaces,
            create_mesh,
            deform_domain!,
            radial_dimension,
            read_gmsh,
            read_cells,
            save_cells,
            read_surfaces,
            save_surfaces,
            read_mesh,
            save_mesh,
            AbstractTimeProfile,
            AbstractGradient,
            integral,
            int_F2,
            intervals,
            echotime,
            isconstant,
            diffusion_time,
            diffusion_time_sta,
            j_integral,
            fibonacci_sphere,
            AbstractCallback,
            initialize!,
            update!,
            finalize!
        """))
end

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
