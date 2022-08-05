"""
    plot_field(femesh, ξ)

Plot field `ξ` on the finite element mesh.
"""
function plot_field end

function plot_field(
    femesh::FEMesh{T,2},
    ξ,
    compartments = 1:length(femesh.points),
) where {T}
    ξ_cmpts = split_field(femesh, ξ)
    fig = Figure()
    ax = Axis(fig[1, 1])
    for icmpt ∈ compartments
        elements = femesh.elements[icmpt]
        points = femesh.points[icmpt]
        colors = abs.(ξ_cmpts[icmpt])
        mesh!(ax, points', elements'; color = colors, shading = false)
    end
    fig
end

function plot_field(
    femesh::FEMesh{T,3},
    ξ,
    compartments = 1:length(femesh.points),
) where {T}
    nboundary = size(femesh.facets, 2)
    ξ_cmpts = split_field(femesh, ξ)
    scene = nothing
    first = true
    for icmpt ∈ compartments, iboundary = 1:nboundary
        facets = femesh.facets[icmpt, iboundary]
        if !isempty(facets)
            points = femesh.points[icmpt]
            colors = abs.(ξ_cmpts[icmpt])
            # colors = fill(1.0 * iboundary, length(colors))
            # colors = iboundary / nboundary
            if first
                scene = mesh(points', facets'; color = colors, shading = false)
                first = false
            else
                mesh!(points', facets'; color = colors, shading = false)
            end
        end
    end
    scene
end
