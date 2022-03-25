"""
    AbstractCallback

Callback to be performed after every time step when solving the BTPDE. The callback is
initialized before time stepping using the `initialize!` function, updated after every step
using `update!`, and finalized after time stepping using `finalize!`.
"""
abstract type AbstractCallback end

"""
    Printer(; nupdate = 1)

Print time stepping information to the console.
"""
Base.@kwdef mutable struct Printer <: AbstractCallback
    verbosity::Int = 2
    nupdate::Int = 1
    n::Int = 1
end

"""
    VTKWriter(; nupdate = 1, dir = "output", filename = "solution")

Write magnetization field to a VTK file after every `nupdate` timestep. The files are stored in a
ParaView data collection file (PVD). The magnetization time series may be visualized in
ParaView by opening the file `"\$dir/\$filename.pvd"`. The compartments are labeled.
"""
Base.@kwdef mutable struct VTKWriter <: AbstractCallback
    nupdate::Int = 1
    n::Int = 1
    ifile::Int = 1
    dir::String = "output/time_stepping"
    filename::String = "solution"
    pvd::WriteVTK.CollectionFile = paraview_collection("")
end

"""
    Plotter(; nupdate = 1)

Plot the evolution of the BTPDE during time stepping. This requires loading the `GLMakie`
plotting backend (`]add GLMakie; using GLMakie`). The plot is updated every `nupdate` time
step. The resulting figure contains a plot of the time profile, total signal attenuation,
and magnetization field (complex magnitude and phase shift).
"""
Base.@kwdef mutable struct Plotter{T} <: AbstractCallback
    nupdate::Int = 1
    n::Int = 1
    t::Observable{Vector{T}} = Observable(T[])
    f::Observable{Vector{T}} = Observable(T[])
    g⃗::Observable{Vector{Vec3f}} = Observable([Vec3f(0, 0, 0)])
    g⃗_hist::Observable{Vector{Vec3f}} = Observable([Vec3f(0, 0, 0)])
    ξ::Observable{Vector{Complex{T}}} = Observable(Complex{T}[])
    magnitude::Matrix{Observable{Vector{T}}} = [Observable(Vector{T}());;]
    phase::Matrix{Observable{Vector{T}}} = [Observable(Vector{T}());;]
    attenuation::Observable{Vector{T}} = Observable(T[])
    S₀::T = one(T)
    fig::Figure = Figure()
end
