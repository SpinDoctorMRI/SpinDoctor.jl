"""
    split_mesh(mesh_all)

Split global mesh into compartments.
"""
function split_mesh(mesh_all)
    # Extract mesh_all in global numbering system (tag: "all")
    points_all = mesh_all.points
    facets_all = mesh_all.facets
    facetmarkers_all = mesh_all.facetmarkers
    elements_all = mesh_all.elements
    elementmarkers_all = mesh_all.elementmarkers

    dim, npoint = size(mesh_all.points)
    ncompartment = maximum(elementmarkers_all)
    nboundary = maximum(facetmarkers_all)
    elements = [elements_all[:, elementmarkers_all.==icmpt] for icmpt = 1:ncompartment]

    # Map global to local points
    point_map = sort.(unique.(elements))
    points = [points_all[:, point_map[icmpt]] for icmpt = 1:ncompartment]

    # Assign each boundary facet to the compartments that have it as an element face
    # (an interface belongs to both of its compartments), and renumber locally
    facets = [zeros(Int, dim, 0) for _ = 1:ncompartment, _ = 1:nboundary]
    point_map_inv = zeros(Int, npoint)
    for icmpt = 1:ncompartment
        point_map_inv[point_map[icmpt]] = 1:length(point_map[icmpt])
        inside = facet_elements(elements[icmpt], facets_all) .> 0
        for iboundary = 1:nboundary
            sel = inside .& (facetmarkers_all .== iboundary)
            if any(sel)
                facets[icmpt, iboundary] = point_map_inv[facets_all[:, sel]]
            end
        end
        elements[icmpt] = point_map_inv[elements[icmpt]]
        point_map_inv[point_map[icmpt]] .= 0
    end

    FEMesh(; point_map, points, facets, elements)
end
