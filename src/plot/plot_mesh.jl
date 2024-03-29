"""
    plot_mesh(femesh, compartments = 1:ncompartment)

Plot finite element mesh, with a subset of the compartments.
"""
function plot_mesh end

function plot_mesh(femesh::FEMesh{T,2}, compartments = 1:length(femesh.points)) where {T}
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y")
    ax.aspect = DataAspect()
    for (i, icmpt) ∈ enumerate(compartments)
        points = femesh.points[icmpt]
        elements = femesh.elements[icmpt]
        poly!(ax, points', elements'; color = Cycled(i), strokewidth = 1, shading = false)
    end
    fig
end

function plot_mesh(
    femesh::FEMesh{T,3},
    compartments = 1:length(femesh.points),
    boundaries = 1:size(femesh.facets, 2),
) where {T}
    scene = Figure()
    first = true
    for icmpt ∈ compartments, iboundary ∈ boundaries
        facets = femesh.facets[icmpt, iboundary]
        if !isempty(facets)
            points = femesh.points[icmpt]
            color = points[3, :]
            if first
                # scene = mesh(points', facets'; color, shading = false)
                scene = poly(points', facets'; color, strokewidth = 1, shading = false)
                first = false
            else
                # mesh!(points', facets'; color, shading = false)
                poly!(points', facets'; color, strokewidth = 1, shading = false)
            end
        end
    end
    scene
end
