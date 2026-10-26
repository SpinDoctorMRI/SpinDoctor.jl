function plot_surfaces(surfaces, boundaries = 1:maximum(surfaces.facetmarkers))
    (; points, facets, facetmarkers) = surfaces
    dim, npoint = size(points)
    if dim == 2
        fig = Figure()
        ax = Axis(fig[1, 1])
        i = 1
        for b in boundaries
            f = facets[:, facetmarkers .== b]
            if !isempty(f)
                segs = [
                    (Point2f(points[:, f[1]]), Point2f(points[:, f[2]])) for f in eachcol(f)
                ]
                linesegments!(ax, segs, color = Cycled(i))
                i += 1
            end
        end
        scene = fig
    elseif dim == 3
        scene = Figure()
        first = true
        for b in boundaries
            f = facets[:, facetmarkers .== b]
            if !isempty(f)
                color = points[3, :]
                if first
                    scene = poly(
                        points',
                        f';
                        color = (:blue, 0.2),
                        strokewidth = 1,
                        shading = NoShading,
                    )
                    first = false
                else
                    poly!(
                        points',
                        f';
                        color = (:blue, 0.2),
                        strokewidth = 1,
                        shading = NoShading,
                    )
                end
            end
        end
    end
    scene
end
