# Functions implemented in package extensions

# Plotting, implemented in the Makie extension, loaded together with a Makie backend (e.g.
# `using CairoMakie` or `using GLMakie`)

"""
    plot_surfaces(surfaces; boundaries = 1:nboundary, title = "", size = (600, 500))

Plot the surface triangulation `surfaces` (as returned by [`create_geometry`](@ref)), with
one color per boundary. Requires a Makie backend.
"""
function plot_surfaces end

"""
    plot_mesh(femesh; title = "", decorations = true, size = (600, 500), kwargs...)
    plot_mesh!(ax, femesh; compartments, colors, alpha = 1, strokewidth, strokecolor)

Plot the finite element mesh `femesh`, with one color per compartment: the elements in 2D,
and the surfaces of the compartments in 3D (each interface once). Only the given
`compartments` are shown. `colors` and `alpha` (transparency) may be given for each
compartment (`alpha` also as one value). With `strokewidth > 0`, the edges of the elements
or facets are drawn.

`plot_mesh` creates a figure (without axis decorations if `decorations = false`), and
`plot_mesh!` plots into the axis `ax` (e.g. an `Axis3` with `aspect = :data`). Requires a
Makie backend.
"""
function plot_mesh end

"""
    plot_mesh!(ax, femesh; compartments, colors, alpha = 1, strokewidth, strokecolor)

Plot the finite element mesh `femesh` into the axis `ax`, with the keywords of
[`plot_mesh`](@ref).
"""
function plot_mesh! end

"""
    plot_field(femesh, u; part = :abs, title = "", decorations = true, size, kwargs...)
    plot_field!(ax, femesh, u; part = :abs, compartments, colormap, colorrange)

Plot the complex field `u` (e.g. the magnetization) on the finite element mesh: its
magnitude (`part = :abs`), real or imaginary part (`:real`, `:imag`) or phase (`:phase`),
on the elements in 2D and on the surfaces of the compartments in 3D. `u` may be an
`Observable`, for animations.

`plot_field` creates a figure with a colorbar, and `plot_field!` plots into the axis `ax`
and returns the plot (e.g. for `Colorbar(position, plot)`). Requires a Makie backend.
"""
function plot_field end

"""
    plot_field!(ax, femesh, u; part = :abs, compartments, colormap, colorrange)

Plot the complex field `u` into the axis `ax`, with the keywords of [`plot_field`](@ref).
Returns the plot (e.g. for `Colorbar(position, plot)`).
"""
function plot_field! end

"""
    plot_hardi(directions, values; title = "", label = "", size, colormap = :viridis)
    plot_hardi!(ax, directions, values; colormap = :viridis, colorrange)

Plot `values` (e.g. signals or ADCs) in the given `directions` (columns of unit vectors) as
a surface (3D) or a closed curve (2D) at the distances `values` from the origin, colored by
the values. Requires a Makie backend.
"""
function plot_hardi end

"""
    plot_hardi!(ax, directions, values; colormap = :viridis, colorrange)

Plot `values` in the given `directions` into the axis `ax` (see [`plot_hardi`](@ref)).
"""
function plot_hardi! end

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
        if f in (
            plot_surfaces,
            plot_mesh,
            plot_mesh!,
            plot_field,
            plot_field!,
            plot_hardi,
            plot_hardi!,
            initialize_plotter!,
            update_plotter!,
        )
            print(
                io,
                "\nPlotting requires a Makie backend. Load one with `using CairoMakie` or `using GLMakie`.",
            )
        elseif f in (save_geometry, load_geometry)
            print(io, "\nSaving and loading geometries requires `using JLD2`.")
        end
    end
end
