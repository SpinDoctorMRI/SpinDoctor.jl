# The plotting functions are implemented in the Makie extension, loaded together with a
# Makie backend (e.g. `using CairoMakie` or `using GLMakie`).

"""
    plot_surfaces(surfaces, boundaries = 1:nboundary)

Plot the surface triangulation `surfaces` (as returned by [`create_geometry`](@ref)), with a
subset of the boundaries. Requires a Makie backend.
"""
function plot_surfaces end

"""
    plot_mesh(femesh, compartments = 1:ncompartment[, boundaries])

Plot finite element mesh, with a subset of the compartments (and boundaries in 3D).
Requires a Makie backend.
"""
function plot_mesh end

"""
    plot_field(femesh, ξ, compartments = 1:ncompartment)

Plot the magnitude of the field `ξ` on the finite element mesh. Requires a Makie backend.
"""
function plot_field end

"""
    plot_hardi(directions, values)

Plot `values` (e.g. signals or ADCs) as radii in the given `directions`. Requires a Makie
backend.
"""
function plot_hardi end

# Implemented in the Makie extension
function initialize_plotter! end
function update_plotter! end

# Saving and loading geometries in JLD2 files (implemented in the JLD2 extension, loaded
# with `using JLD2`)

"""
    save_geometry(file, setup, mesh, surfaces, cells)

Save a geometry created by [`create_geometry`](@ref) to the JLD2 file `file`, together with
a description of `setup`. Only plain arrays are stored, so the file does not depend on the
SpinDoctor version. Requires `using JLD2`.
"""
function save_geometry end

"""
    (; setup, mesh, surfaces, cells) = load_geometry(file)

Load a geometry saved with [`save_geometry`](@ref). The `setup` is returned as the
description string that was saved. Requires `using JLD2`.
"""
function load_geometry end

# Hints for functions implemented in package extensions
function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        f = exc.f
        if f ∈ (
            plot_surfaces,
            plot_mesh,
            plot_field,
            plot_hardi,
            initialize_plotter!,
            update_plotter!,
        )
            print(
                io,
                "\nPlotting requires a Makie backend. Load one with `using CairoMakie` or `using GLMakie`.",
            )
        elseif f ∈ (save_geometry, load_geometry)
            print(io, "\nSaving and loading geometries requires `using JLD2`.")
        end
    end
end
