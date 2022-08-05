"""
    read_mesh(filename)
    
Read mesh from Tetgen format file.
"""
function read_mesh(filename)
    @info "Reading Tetgen FE mesh from " * filename

    # Read points
    dim, points = open(filename * ".node", "r") do io
        npoint, dim, = parse.(Int, split(readline(io)))
        points = zeros(dim, npoint)
        for ipoint = 1:npoint
            points[:, ipoint] = parse.(Float64, split(readline(io))[2:4])
        end
        points
    end

    # Read facets and their associated boundaries
    facets, facetmarkers = open(filename * ".face", "r") do io
        nfacet = parse(Int, split(readline(io))[1])
        facets = zeros(Int, dim, nfacet)
        facetmarkers = zeros(Int, nfacet)
        for ifacet = 1:nfacet
            tmp = parse.(Int, split(readline(io)))
            facets[:, ifacet] = tmp[2:1+dim]
            facetmarkers[ifacet] = tmp[2+dim]
        end
        facets, facetmarkers
    end

    # Read elements and their associated compartments
    elements, elementmarkers = open(filename * ".ele", "r") do io
        nelement = parse(Int, split(readline(io))[1])
        elements = zeros(Int, dim + 1, nelement)
        elementmarkers = zeros(Int, nelement)
        for ielement = 1:nelement
            tmp = parse.(Int, split(readline(io)))
            elements[:, ielement] = tmp[2:2+dim]
            elementmarkers[ielement] = tmp[3+dim]
        end
        elements, elementmarkers
    end

    (; points, facets, facetmarkers, elements, elementmarkers)
end
