using SpinDoctor
using LinearAlgebra
using Test


@time @testset "Create cells" begin
    include("test_create_cells.jl")
end

@time @testset "Create surface triangulation" begin
    include("test_create_surface_triangulation.jl")
end
