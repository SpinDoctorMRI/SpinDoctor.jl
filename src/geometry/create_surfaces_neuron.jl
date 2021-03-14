"""
    create_surfaces_neuron(filename, setup)

Create neuron surface mesh.
A neuron surface mesh is loaded or create and loaded. An ECS can be added.
"""
function create_surfaces_neuron(filename, setup::Setup)

    # Number of points to discretize space for creating tight wrap ECS
    ndiscretize = 100

    ecs_shape = setup.geometry[:ecs_shape]
    ecs_ratio = setup.geometry[:ecs_ratio]

    if !isfile(filename * "_elements.txt") || !isfile(filename * "_nodes.txt")
        gmsh_to_fem(filename)
    end
    println("Reading from neuron FE mesh from " * filename)

    # Read points
    nline = open(filename * "_nodes.txt", "r") do io
        countlines(io)
    end
    points = open(filename * "_nodes.txt", "r") do io
        points = zeros(3, nline)
        for iline = 1:nline
            points[:, iline] = parse.(Float64, split(readline(io)))
        end
        points
    end

    # Read elements
    nline = open(filename * "_elements.txt", "r") do io
        countlines(io)
    end
    elements = open(filename * "_elements.txt", "r") do io
        elements = zeros(Int, 4, nline)
        for iline = 1:nline
            elements[:, iline] = parse.(Int, split(readline(io)))
        end
        elements
    end

    # Extract boundary facets from elements (remove interior facets)
    facets = [elements[[1, 2, 3], :] elements[[1, 2, 4], :] elements[[1, 3, 4], :] elements[
        [2, 3, 4],
        :,
    ]]
    sort!(facets, dims = 1)
    facets_unq = unique(facets, dims = 2)
    ikeep = [count(f -> f == fu, eachcol(facets)) == 1 for fu ∈ eachcol(facets_unq)]
    facets = facets_unq[:, ikeep]

    volumes, centers = get_mesh_volumes(points, elements)
    centermass = centers * volumes / sum(volumes)

    facetmarkers = ones(Int, size(facets, 2))
    regions = centermass

    pmin = minimum(points, dims = 2)
    pmax = maximum(points, dims = 2)

    if ecs_shape != "no_ecs"

        # ecs_gap = ecs_ratio * min(pmax - pmin)
        ecs_gap = ecs_ratio * 10

        if ecs_shape == "box"
            points_ecs = [
                pmin(1) pmax(1) pmax(1) pmin(1) pmin(1) pmax(1) pmax(1) pmin(1)
                pmin(2) pmin(2) pmax(2) pmax(2) pmin(2) pmin(2) pmax(2) pmax(2)
                pmin(3) pmin(3) pmin(3) pmin(3) pmax(3) pmax(3) pmax(3) pmax(3)
            ]
            facets_ecs = [
                1 1 1 1 2 2 3 3 4 4 5 5
                2 3 2 6 3 7 4 8 1 5 6 7
                3 4 6 5 7 6 8 7 5 8 7 8
            ]
        elseif ecs_shape == "convex_hull"
            error("Not implemented")
        elseif ecs_shape == "tight_wrap"
            error("Not implemented")
        end

        facetmarkers_ecs = 2 * ones(Int, size(facets_ecs, 2))
        _, ind = sort(points[1, :])
        ind = ind[1]
        regions_ecs = points[:, ind] - ecs_gap / 10 * [1; 0; 0]

        npoint_out = size(points, 2)
        points = [points points_ecs]
        facets = [facets facets_ecs + npoint_out]
        append!(facetmarkers, facetmarkers_ecs)
        append!(regions, regions_ecs)
    end


    # Return named tuple
    (; points, facets, facetmarkers, regions)
end
