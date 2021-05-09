function split_mesh(mesh_all)

    # Extract mesh_all in global numbering system (tag: "all")
    points_all = mesh_all.points
    facets_all = mesh_all.facets
    facetmarkers_all = mesh_all.facetmarkers
    elements_all = mesh_all.elements
    elementmarkers_all = mesh_all.elementmarkers

    ncompartment = maximum(elementmarkers_all)
    nboundary = maximum(facetmarkers_all)
    npoint_all = size(points_all, 2)
    elements = [elements_all[:, elementmarkers_all.==icmpt] for icmpt = 1:ncompartment]
    point_map = sort.(unique.(elements))
    ncompartment_inds = length.(point_map)

    boundary_facets =
        [facets_all[:, facetmarkers_all.==iboundary] for iboundary = 1:nboundary]
    boundary_inds = sort.(unique.(boundary_facets))
    nboundary_inds = length.(boundary_inds)

    points = [points_all[:, point_map[icmpt]] for icmpt = 1:ncompartment]

    facets = fill(zeros(Int, 3, 0), ncompartment, nboundary)
    for icmpt = 1:ncompartment
        println("Separating compartment $icmpt of $ncompartment")
        # old_new = Pair.(point_map[icmpt], 1:length(point_map[icmpt]))

        # replace!(elements[icmpt], old_new...)

        replace!(x -> findfirst(x .== point_map[icmpt]), elements[icmpt])

        # for ipoint = 1:ncompartment_inds[icmpt]
        #     elements[icmpt][elements[icmpt].==point_map[icmpt][ipoint]] .= ipoint
        # end

        for iboundary = 1:nboundary
            println("Separating boundary $iboundary of $nboundary")
            if all(boundary_facets[iboundary] .∈ [point_map[icmpt]])
                # facets[icmpt, iboundary] = replace(boundary_facets[iboundary], old_new...)

                facets[icmpt, iboundary] = replace(
                    x -> findfirst(x .== point_map[icmpt]),
                    boundary_facets[iboundary],
                )

                # facets[icmpt, iboundary] = deepcopy(boundary_facets[iboundary])
                # for ipoint = 1:ncompartment_inds[icmpt]
                #     facets[icmpt, iboundary][facets[
                #         icmpt,
                #         iboundary,
                #     ].==point_map[icmpt][ipoint]] .= ipoint
                # end
            end
        end
    end

    FEMesh(; point_map, points, facets, elements)
end
