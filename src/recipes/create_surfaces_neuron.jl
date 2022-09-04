function create_surfaces(setup::NeuronSetup{T}, filename) where {T}
    (; ecs) = setup

    if !isfile(filename * "_elements.txt") || !isfile(filename * "_nodes.txt")
        gmesh2fem(filename)
    end
    @info "Reading neuron FE mesh from " * filename

    # Read points
    nline = open(filename * "_nodes.txt", "r") do io
        countlines(io)
    end
    points = open(filename * "_nodes.txt", "r") do io
        points = zeros(3, nline)
        for iline = 1:nline
            points[:, iline] = parse.(T, split(readline(io)))
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
    # FIXME: The below statement is slow
    ikeep = [count(f -> f == fu, eachcol(facets)) == 1 for fu ∈ eachcol(facets_unq)]
    facets = facets_unq[:, ikeep]

    volumes, centers = get_mesh_volumes(points, elements)
    centermass = centers * volumes / sum(volumes)

    facetmarkers = ones(Int, size(facets, 2))
    regions = centermass

    # Create ECS around points
    ecs_surface = create_surfaces(ecs, points)

    npoint = size(points, 2)
    points = [points ecs_surface.points]
    facets = [facets ecs_surface.facets .+ npoint]
    append!(facetmarkers, fill(2, size(ecs_surface.facets, 2)))
    regions = [regions ecs_surface.region]

    (; points, facets, facetmarkers, regions)
end
