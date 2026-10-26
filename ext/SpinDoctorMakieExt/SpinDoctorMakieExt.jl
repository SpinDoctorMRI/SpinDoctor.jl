"""
Plotting functions and the [`Plotter`](@ref) callback, available when a Makie backend
(e.g. `CairoMakie` or `GLMakie`) is loaded.
"""
module SpinDoctorMakieExt

using LinearAlgebra
using Makie
using SpinDoctor
using SpinDoctor: FEMesh, Plotter, convexhull, echotime
import SpinDoctor: plot_surfaces, plot_mesh, plot_field, plot_hardi
import SpinDoctor: plot_mesh!, plot_field!, plot_hardi!
import SpinDoctor: initialize_plotter!, update_plotter!

include("meshes.jl")
include("plot_hardi.jl")
include("plot_surfaces.jl")
include("plotter.jl")

end
