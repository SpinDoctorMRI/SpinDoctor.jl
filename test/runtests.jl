# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/SpinDoctor.jl")
    using .SpinDoctor
end

using SpinDoctor
using LinearAlgebra
using QuadGK
using Test
using Aqua

include("testutils.jl")

# Run through typical workflow
@testset "Workflow" begin
    include("workflow.jl")
end

@testset "Aqua" begin
    Aqua.test_all(
        SpinDoctor;
        ambiguities = false,
        project_toml_formatting = false, # https://github.com/JuliaTesting/Aqua.jl/issues/72
    )
end
