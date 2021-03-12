function split_mesh(mesh_all, setup)
    @unpack boundary_markers, compartments, boundaries = setup.pde
    nboundary, ncompartment = size(boundary_markers)

    # Extract mesh_all in global numbering system (tag: "all")
    points_all = mesh_all.points
    facets_all = mesh_all.facets
    facetmarkers_all = mesh_all.facetmarkers
    elements_all = mesh_all.elements
    elementmarkers_all = mesh_all.elementmarkers

    ncompartment = maximum(elementmarkers_all)
    nboundary = maximum(facetmarkers_all)
    npoint_all = size(points_all, 2)

    elements = [elements_all[:, elementmarkers_all .== icmpt] for icmpt = 1:ncompartment]
    cmpt_inds = sort.(unique.(elements))
    ncompartment_inds = length.(cmpt_inds)

    facets_boundary = [facets_all[:, facetmarkers_all .== iboundary] for iboundary = 1:nboundary]
    boundary_inds = sort.(unique.(facets_boundary))
    nboundary_inds = length.(boundary_inds)

    points = [points_all[:, cmpt_inds[icmpt]] for icmpt = 1:ncompartment]

    facets = fill(zeros(Int, 3, 0), ncompartment, nboundary)
    for icmpt = 1:ncompartment
        for ipoint = 1:ncompartment_inds[icmpt]
            elements[icmpt][elements[icmpt] .== cmpt_inds[icmpt][ipoint]] .= ipoint
        end
        for iboundary = 1:nboundary
            if boundary_markers[icmpt, iboundary]
                facets[icmpt, iboundary] = deepcopy(facets_boundary[iboundary])
                for ipoint = 1:ncompartment_inds[icmpt]
                    facets[icmpt, iboundary][facets[icmpt, iboundary] .== cmpt_inds[icmpt][ipoint]] .= ipoint
                end
            end
        end
    end

    (; ncompartment, nboundary, cmpt_inds, points, facets, elements, boundary_markers)
end
