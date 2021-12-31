"""
    plot_mesh(femesh)

Plot finite element mesh.
"""
function plot_mesh(femesh::FEMesh)
    ncompartment, nboundary = size(femesh.facets)
    scene = nothing
    first = true
    for icmpt = 1:ncompartment, iboundary = 1:nboundary
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
