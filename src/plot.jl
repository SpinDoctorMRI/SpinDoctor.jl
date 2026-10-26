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
        end
    end
end
