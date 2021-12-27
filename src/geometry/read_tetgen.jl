"""
    read_tetgen
    
Read mesh from Tetgen.
"""
function read_tetgen(filename)

    @info "Reading Tetgen FE mesh from " * filename

    # Read points
    points = open(filename * ".node", "r") do io
        npoint = parse(Int, split(readline(io))[1])
        points = zeros(3, npoint)
        for ipoint = 1:npoint
            points[:, ipoint] = parse.(Float64, split(readline(io))[2:4])
        end
        points
    end

    # Read facets and their associated boundaries
    facets, facetmarkers = open(filename * ".face", "r") do io
        nfacet = parse(Int, split(readline(io))[1])
        facets = zeros(Int, 3, nfacet)
        facetmarkers = zeros(Int, nfacet)
        for ifacet = 1:nfacet
            tmp = parse.(Int, split(readline(io)))
            facets[:, ifacet] = tmp[2:4]
            facetmarkers[ifacet] = tmp[5]
        end
        facets, facetmarkers
    end

    # Read elements and their associated compartments
    elements, elementmarkers = open(filename * ".ele", "r") do io
        nelement = parse(Int, split(readline(io))[1])
        elements = zeros(Int, 4, nelement)
        elementmarkers = zeros(Int, nelement)
        for ielement = 1:nelement
            tmp = parse.(Int, split(readline(io)))
            elements[:, ielement] = tmp[2:5]
            elementmarkers[ielement] = tmp[6]
        end
        elements, elementmarkers
    end

    (; points, facets, facetmarkers, elements, elementmarkers)
end
