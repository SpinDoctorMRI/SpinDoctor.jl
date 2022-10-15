"""
    read_surfaces(filename)

Read surface triangulation from file.
"""
function read_surfaces(filename)
    @info "Loading surface triangulation from " * filename

    # Get list of nodes from .node file
    points, dim = open(filename * ".node", "r") do io
        readline(io)
        readline(io)
        npoint, dim = parse.(Int, split(readline(io))[1:2])
        readline(io)
        points = zeros(dim, npoint)
        for ipoint = 1:npoint
            tmp = parse.(Float64, split(readline(io)))
            points[:, ipoint] = tmp[2:2+dim-1]
        end
        points, dim
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
        facets = zeros(Int, dim, nfacet)
        facetmarkers = zeros(Int, nfacet)
        for ifacet = 1:nfacet
            tmp = parse.(Int, split(readline(io)))
            facetmarkers[ifacet] = tmp[3]
            tmp = parse.(Int, split(readline(io)))
            facets[:, ifacet] = tmp[2:2+dim-1]
        end

        # Read list of holes (empty)
        readline(io)
        readline(io)

        # Read list of interior points with their corresponding compartment
        readline(io)
        nregion = parse(Int, split(readline(io))[1])
        regions = zeros(dim, nregion)
        for iregion = 1:nregion
            tmp = parse.(Float64, split(readline(io)))
            regions[:, iregion] = tmp[2:2+dim-1]
        end

        facets, facetmarkers, regions
    end

    (; points, facets, facetmarkers, regions)
end
