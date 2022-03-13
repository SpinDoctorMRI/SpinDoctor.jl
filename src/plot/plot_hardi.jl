function plot_hardi(directions, values)
    facets, points = convexhull(directions)
    for (i, v) ∈ enumerate(values)
        points[:, i] .*= v
    end
    mesh(points', facets', color = values, shading = false)
end
