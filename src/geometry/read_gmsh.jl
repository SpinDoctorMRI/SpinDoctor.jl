"""
    points, elements = read_gmsh(filename)

Read the points and tetrahedra of a volume mesh in the Gmsh format (version 2, as in the
RealNeuronMeshes repository).
"""
function read_gmsh(filename)
    open(filename, "r") do io
        # Header
        for _ = 1:4
            readline(io)
        end

        # Points
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
            elements[:, ielement] = parse.(Int, split(readline(io))[6:9])
        end

        points, elements
    end
end
