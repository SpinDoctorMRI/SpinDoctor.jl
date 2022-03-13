# # High angular resolution diffusion imaging
#
# In this example we will compute the signal using the same gradient sequence in many
# different directions.

using SpinDoctor
using LinearAlgebra

if haskey(ENV, "GITHUB_ACTIONS")
    using CairoMakie
else
    using GLMakie
end

# LSP indexing solution                                                          #src
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983 #src
if isdefined(@__MODULE__, :LanguageServer)                                       #src
    include("../src/SpinDoctor.jl")                                              #src
    using .SpinDoctor                                                            #src
end                                                                              #src

# TODO: Write example
