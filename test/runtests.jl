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
using Random
using JLD2

include("testutils.jl")

@testset "SpinDoctor" begin
    @testset "Gradient sequences" begin
        include("datatypes/sequences.jl")
    end

    @testset "Checks" begin
        include("checks.jl")
    end

    # Run through typical workflow
    @testset "Workflow" begin
        include("workflow.jl")
    end

    @testset "Aqua" begin
        Aqua.test_all(SpinDoctor; ambiguities = false)
    end
end
