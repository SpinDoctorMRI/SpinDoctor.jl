"""
    gmesh2fem(name)

Extract points and elements from Gmesh file `name.msh`. This will create the two files
`name_nodes.txt` and `name_elements.txt` in the same directory.
"""
function gmesh2fem(name)
    # Read points and elements
    mshname = name * ".msh"
    @info "Reading nodes and elements from " * mshname
    points, elements = open(mshname, "r") do io
        # Points
        for i = 1:4
            readline(io)
        end
        npoint = parse(Int, readline(io))
        points = zeros(3, npoint)
        for ipoint = 1:npoint
            points[:, ipoint] = parse.(Float64, split(readline(io))[2:4])
        end
        readline(io)

        # Elements
        readline(io)
        nelement = parse(Int, readline(io))
        elements = zeros(Int, 4, nelement)
        for ielement = 1:nelement
            tmp = parse.(Int, split(readline(io)))
            elements[:, ielement] = tmp[6:9]
        end

        points, elements
    end

    # Write nodes
    nodename = name * "_nodes.txt"
    @info "Writing nodes to " * nodename
    npoint = size(points, 2)
    open(nodename, "w") do io
        for i = 1:npoint
            for j = 1:3
                write(io, "  ")
                write(io, string(points[j, i]))
            end
            write(io, "\n")
        end
    end

    # Write elements
    elementname = name * "_elements.txt"
    @info "Writing elements to " * elementname
    nelement = size(elements, 2)
    open(elementname, "w") do io
        for i = 1:nelement
            for j = 1:4
                write(io, "  ")
                write(io, string(elements[j, i]))
            end
            write(io, "\n")
        end
    end

    nothing
end
