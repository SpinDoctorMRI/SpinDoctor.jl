function plot_surfaces(
    surfaces;
    boundaries = 1:maximum(surfaces.facetmarkers),
    title = "",
    size = (600, 500),
)
    (; points) = surfaces
    facets, facetmarkers = SpinDoctor.triangulate_facets(surfaces)
    dim = Base.size(points, 1)
    fig = Figure(; size)
    ax = mesh_axis(fig[1, 1], dim; title)
    for (i, b) in enumerate(boundaries)
        f = facets[:, facetmarkers .== b]
        isempty(f) && continue
        color = Makie.wong_colors()[mod1(i, 7)]
        if dim == 2
            segments =
                [(Point2f(points[:, e[1]]), Point2f(points[:, e[2]])) for e in eachcol(f)]
            linesegments!(ax, segments; color)
        else
            poly!(
                ax,
                points',
                f';
                color = (color, 0.3),
                strokewidth = 0.5,
                strokecolor = color,
            )
        end
    end
    fig
end
