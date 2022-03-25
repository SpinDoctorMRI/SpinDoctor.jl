"""
    plot_mesh(femesh, compartments = 1:ncompartment)

Plot finite element mesh, with a subset of the compartments.
"""
function plot_mesh(femesh::FEMesh, compartments = 1:length(femesh.points))
    nboundary = size(femesh.facets, 2)
    scene = nothing
    first = true
    for icmpt ∈ compartments, iboundary = 1:nboundary
        facets = femesh.facets[icmpt, iboundary]
        points = femesh.points[icmpt]
        color = points[3, :]
        # color = fill(1.0 * iboundary, length(colors))
        # color = iboundary / nboundary
        if first
            # scene = mesh(points', facets'; color, shading = false)
            scene = poly(points', facets'; color, strokewidth = 1, shading = false)
            first = false
        else
            # mesh!(points', facets'; color, shading = false)
            poly!(points', facets'; color, strokewidth = 1, shading = false)
        end
    end
    scene
end
