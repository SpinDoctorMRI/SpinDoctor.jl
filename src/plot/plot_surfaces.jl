"""
    plot_surfaces(femesh, compartments = 1:ncompartment)

Plot surfaces.
"""
function plot_surfaces(surfaces, boundaries = 1:maximum(surfaces.facetmarkers))
    (; points, facets, facetmarkers) = surfaces
    dim = size(points, 1)
    if dim == 2
        fig = Figure()
        ax = Axis(fig[1, 1])
        for (i, b) ∈ enumerate(boundaries)
            f = facets[:, facetmarkers .== b]
            segs = [(Point2f(points[:, f[1]]), Point2f(points[:, f[2]])) for f ∈ eachcol(f)]
            linesegments!(ax, segs, color = Cycled(i))
        end
        scene = fig
    elseif dim == 3
        scene = nothing
        first = true
        for b ∈ boundaries
            f = facets[:, facetmarkers .== b]
            if !isempty(f)
                color = points[3, :]
                if first
                    scene = poly(points', f'; color, strokewidth = 1, shading = false)
                    first = false
                else
                    poly!(points', f'; color, strokewidth = 1, shading = false)
                end
            end
        end
    end
    scene
end
