using SpinDoctor, Test, SafeTestsets

@time @safetestset "Create cells" begin
    include("test_create_cells.jl")
end
@time @safetestset "Create surface triangulation" begin
    include("test_create_surface_triangulation.jl")
end
