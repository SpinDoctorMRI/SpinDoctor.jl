function split_mesh(domain, mesh)
    @unpack boundary_markers, compartments, boundaries, ncmpt, nboundary = domain

    # Extract mesh in global numbering system (tag: "all")
    points_all = mesh.pointlist
    facets_all = Int.(mesh.trifacelist)
    elements_all = Int.(mesh.tetrahedronlist)
    cmptmarkers_all = Int.(mesh.tetrahedronattributelist[:])
    boundarymarkers_all = Int.(mesh.trifacemarkerlist)

    ncmpt = maximum(cmptmarkers_all)
    nboundary = maximum(boundarymarkers_all)
    npoint_all = size(points_all, 2)

    elements = [elements_all[:, cmptmarkers_all .== icmpt] for icmpt = 1:ncmpt]
    cmpt_inds = sort.(unique.(elements))
    ncmpt_inds = length.(cmpt_inds)

    facets_boundary = [facets_all[:, boundarymarkers_all .== iboundary] for iboundary = 1:nboundary]
    boundary_inds = sort.(unique.(facets_boundary))
    nboundary_inds = length.(boundary_inds)

    points = [points_all[:, cmpt_inds[icmpt]] for icmpt = 1:ncmpt]

    facets = fill(zeros(Int, 3, 0), ncmpt, nboundary)
    for icmpt = 1:ncmpt
        for ipoint = 1:ncmpt_inds[icmpt]
            elements[icmpt][elements[icmpt] .== cmpt_inds[icmpt][ipoint]] .= ipoint
        end
        for iboundary = 1:nboundary
            if boundary_markers[icmpt, iboundary]
                facets[icmpt, iboundary] = deepcopy(facets_boundary[iboundary])
                for ipoint = 1:ncmpt_inds[icmpt]
                    facets[icmpt, iboundary][facets[icmpt, iboundary] .== cmpt_inds[icmpt][ipoint]] .= ipoint
                end
            end
        end
    end

    (; ncmpt, nboundary, cmpt_inds, points, facets, elements, boundary_markers)
end
