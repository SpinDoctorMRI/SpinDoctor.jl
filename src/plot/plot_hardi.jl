function plot_hardi(directions, values)
    facets, points = convexhull(directions)
    for (i, v) ∈ enumerate(values)
        points[:, i] .*= v
    end
    poly(points', facets'; strokewidth = 1, color = values, shading = false)
end
