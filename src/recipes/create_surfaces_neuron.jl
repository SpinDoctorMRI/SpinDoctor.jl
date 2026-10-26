function create_surfaces(setup::NeuronSetup, meshfile)
    (; ecs) = setup

    @info "Reading neuron FE mesh from " * meshfile
    points, elements = read_gmsh(meshfile)

    # Boundary facets: the element faces that belong to one element only
    count = Dict{SVector{3,Int},Int}()
    for e ∈ axes(elements, 2), skip = 1:4
        f = sort(SVector{3}(ntuple(k -> elements[k < skip ? k : k + 1, e], Val(3))))
        count[f] = get(count, f, 0) + 1
    end
    facets = reduce(hcat, [Vector(f) for (f, c) ∈ count if c == 1])

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
