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
