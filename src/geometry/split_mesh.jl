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

    npoint = size(mesh_all.points, 2)
    ncompartment = maximum(elementmarkers_all)
    nboundary = maximum(facetmarkers_all)
    elements = [elements_all[:, elementmarkers_all.==icmpt] for icmpt = 1:ncompartment]

    # Map global to local points
    point_map = sort.(unique.(elements))

    # Inverse point map
    point_map_inv = [zeros(Int, npoint) for _ = 1:ncompartment]
    for icmpt = 1:ncompartment
        for ipoint = 1:length(point_map[icmpt])
            point_map_inv[icmpt][point_map[icmpt][ipoint]] = ipoint
        end
    end

    boundary_facets =
        [facets_all[:, facetmarkers_all.==iboundary] for iboundary = 1:nboundary]
    points = [points_all[:, point_map[icmpt]] for icmpt = 1:ncompartment]

    facets = [zeros(Int, 3, 0) for _ = 1:ncompartment, _ = 1:nboundary]
    for icmpt = 1:ncompartment
        # replace!(ipoint -> point_map_inv[icmpt][ipoint], elements[icmpt])
        elements[icmpt] .= point_map_inv[icmpt][elements[icmpt]]

        for iboundary = 1:nboundary
            # f = replace(ipoint -> point_map_inv[icmpt][ipoint], boundary_facets[iboundary])
            f = point_map_inv[icmpt][boundary_facets[iboundary]]
            if all(f .!= 0)
                facets[icmpt, iboundary] = f
            end
        end
    end

    FEMesh(; point_map, points, facets, elements)
end
