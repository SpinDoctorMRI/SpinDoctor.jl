using LinearAlgebra
using QuadGK
using Test

# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/SpinDoctor.jl")
    using .SpinDoctor
else
    using SpinDoctor
end

include("testutils.jl")

# Run through typical workflow
@testset "Workflow" begin
    include("workflow.jl")
end

#= TODO: Restructure tests
# Datatypes
@testset "Setups" begin
    include("datatypes/sequences.jl")
    include("datatypes/experiment.jl")
    include("datatypes/femesh.jl")
    include("datatypes/gradient.jl")
    include("datatypes/model.jl")
    include("datatypes/create_model.jl")
    include("datatypes/setup.jl")
    include("datatypes/get_coefficients.jl")
end

# Geometry
@testset "Geometry" begin
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
end

# Matrix assembly
@testset "Matrix assembly" begin
    include("matrix_assembly/assemble_mass_matrix.jl")
    include("matrix_assembly/assemble_stiffness_matrix.jl")
    include("matrix_assembly/assemble_flux_matrices.jl")
    include("matrix_assembly/assemble_flux_matrix.jl")
    include("matrix_assembly/couple_flux_matrix.jl")
end

# Matrix formalism
@testset "Matrix formalism" begin
    include("matrix_formalism/eig2length.jl")
    include("matrix_formalism/compute_laplace_eig.jl")
    include("matrix_formalism/solve_mf.jl")
end

# Analytical
@testset "Analytical" begin
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
end

# BTPDE
@testset "BTPDE" begin
    include("btpde/solve_btpde.jl")
    include("btpde/solve_btpde_midpoint.jl")
    include("btpde/solve_karger.jl")
end

# ADC
@testset "ADC" begin
    include("adc/solve_hadc.jl")
end

# Postprocessing
@testset "Postprocess" begin
    include("postprocess/fit_adc.jl")
    include("postprocess/fit_tensors.jl")
    include("postprocess/savefield.jl")
end

# Plot
@testset "Plot" begin
    include("plot/plot_mesh.jl")
end
=#
