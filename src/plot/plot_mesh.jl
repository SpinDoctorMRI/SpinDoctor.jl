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
        colors = points[3, :]
        # colors = iboundary / nboundary
        if first
            scene = mesh(points', facets', color = colors, shading = false)
            first = false
        else
            mesh!(points', facets', color = colors, shading = false)
        end
    end
    scene
end
