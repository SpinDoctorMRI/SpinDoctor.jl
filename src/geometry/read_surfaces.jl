"""
    read_surfaces(filename)

Read surface triangulation from file.
"""
function read_surfaces(filename)

    println("Loading surface triangulation from " * filename)

    # Get list of nodes from .node file
    points = open(filename * ".node", "r") do io
        readline(io)
        readline(io)
        npoint = parse(Int, split(readline(io))[1])
        readline(io)
        points = zeros(3, npoint)
        for ipoint = 1:npoint
            vec = parse.(Float64, split(readline(io)))
            points[:, ipoint] = vec[2:4]
        end
        points
    end

    # Read .poly file
    facets, facetmarkers, regions = open(filename * ".poly", "r") do io

        # Read list of holes (refer to separate file)
        readline(io)
        readline(io)
        readline(io)

        # Read list of facets
        readline(io)
        readline(io)
        nfacet = parse(Int, split(readline(io))[1])
        readline(io)
        facets = zeros(Int, 3, nfacet)
        facetmarkers = zeros(Int, nfacet)
        for ifacet = 1:nfacet
            vec = parse.(Int, split(readline(io)))
            facetmarkers[ifacet] = vec[3]
            vec = parse.(Int, split(readline(io)))
            facets[:, ifacet] = vec[2:4]
        end

        # Read list of holes (empty)
        readline(io)
        readline(io)

        # Read list of interior points with their corresponding compartment
        readline(io)
        nregion = parse(Int, split(readline(io))[1])
        regions = zeros(3, nregion)
        for iregion = 1:nregion
            vec = parse.(Float64, split(readline(io)))
            regions[:, iregion] = vec[2:4]
        end

        facets, facetmarkers, regions
    end

    (; points, facets, facetmarkers, regions)
end
