"""
    plot_surfaces(femesh, compartments = 1:ncompartment)

Plot surfaces.
"""
function plot_surfaces(surfaces, boundaries = 1:maximum(surfaces.facetmarkers))
    (; points, facets, facetmarkers) = surfaces
    scene = nothing
    first = true
    for b ∈ boundaries
        f = facets[:, facetmarkers .== b]
        color = points[3, :]
        if first
            scene = poly(points', f'; color, strokewidth = 1, shading = false)
            first = false
        else
            poly!(points', f'; color, strokewidth = 1, shading = false)
        end
    end
    scene
end
