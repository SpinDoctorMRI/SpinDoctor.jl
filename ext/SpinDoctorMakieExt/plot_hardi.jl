# The faces of the convex hull of the directions, in the numbering of the directions
function hull_faces(directions)
    faces, points = convexhull(directions)
    index = map(eachcol(points)) do p
        argmin(j -> sum(abs2, directions[:, j] - p), axes(directions, 2))
    end
    index[faces]
end

function plot_hardi!(ax, directions, values; colormap = :viridis, colorrange = nothing)
    dim = size(directions, 1)
    points = directions .* values'
    if dim == 2
        # Closed curve through the points, in the order of the angles
        order = sortperm(atan.(directions[2, :], directions[1, :]))
        curve = Point2f.(eachcol(points[:, [order; order[1]]]))
        lines!(ax, curve; color = :black)
        scatter!(ax, curve[1:(end-1)]; color = values[order], colormap)
    else
        isnothing(colorrange) && (colorrange = extrema(values))
        mesh!(
            ax,
            points',
            hull_faces(directions)';
            color = values,
            colormap,
            colorrange,
            shading = NoShading,
        )
    end
end

function plot_hardi(
    directions,
    values;
    title = "",
    label = "",
    decorations = true,
    size = (600, 500),
    kwargs...,
)
    dim = Base.size(directions, 1)
    fig = Figure(; size)
    ax = mesh_axis(fig[1, 1], dim; title, decorations)
    p = plot_hardi!(ax, directions, values; kwargs...)
    dim == 3 && Colorbar(fig[1, 2], p; label)
    fig
end
