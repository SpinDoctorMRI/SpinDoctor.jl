function plot_field(
    femesh::FEMesh{T,2},
    u,
    compartments = 1:length(femesh.points),
) where {T}
    u_cmpts = split_field(femesh, u)
    fig = Figure()
    ax = Axis(fig[1, 1])
    for icmpt in compartments
        elements = femesh.elements[icmpt]
        points = femesh.points[icmpt]
        colors = abs.(u_cmpts[icmpt])
        mesh!(ax, points', elements'; color = colors, shading = NoShading)
    end
    fig
end

function plot_field(
    femesh::FEMesh{T,3},
    u,
    compartments = 1:length(femesh.points),
) where {T}
    nboundary = size(femesh.facets, 2)
    u_cmpts = split_field(femesh, u)
    scene = nothing
    first = true
    for icmpt in compartments, iboundary = 1:nboundary
        facets = femesh.facets[icmpt, iboundary]
        if !isempty(facets)
            points = femesh.points[icmpt]
            colors = abs.(u_cmpts[icmpt])
            # colors = fill(1.0 * iboundary, length(colors))
            # colors = iboundary / nboundary
            if first
                scene = mesh(points', facets'; color = colors, shading = NoShading)
                first = false
            else
                mesh!(points', facets'; color = colors, shading = NoShading)
            end
        end
    end
    scene
end
