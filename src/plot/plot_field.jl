"""
    plot_field(femesh, ξ)

Plot field `ξ` on the finite element mesh.
"""
function plot_field(femesh::FEMesh, ξ)
    ncompartment, nboundary = size(femesh.facets)
    ξ_cmpts = split_field(femesh, ξ)
    scene = nothing
    first = true
    for icmpt = 1:ncompartment, iboundary = 1:nboundary
        facets = femesh.facets[icmpt, iboundary]
        points = femesh.points[icmpt]
        colors = abs.(ξ_cmpts[icmpt])
        # colors = fill(1.0 * iboundary, length(colors))
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
