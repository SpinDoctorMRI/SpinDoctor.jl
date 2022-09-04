function plot_hardi(directions, values)
    dim = size(directions, 1)
    if dim == 2
        facets, points = convexhull(directions)
        for (i, v) ∈ enumerate(values)
            points[:, i] .*= v
        end
        segs =
            [(Point2f(points[:, f[1]]), Point2f(points[:, f[2]])) for f ∈ eachcol(facets)]
        fig = Figure()
        ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y")
        linesegments!(ax, segs)
        scatter!(ax, points)
        scatter!(ax, Point2f(0, 0))
        fig
    elseif dim == 3
        facets, points = convexhull(directions)
        for (i, v) ∈ enumerate(values)
            points[:, i] .*= v
        end
        poly(points', facets'; strokewidth = 1, color = values, shading = false)
    end
end
