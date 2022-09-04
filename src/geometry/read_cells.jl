"""
    read_cells(cellfilename)

Read cell configuration (centers and radii).
"""
function read_cells(cellfilename)

    # Read geometry from file
    @info "Loading cell configuration from " * cellfilename

    open(cellfilename, "r") do io
        readline(io)
        ncell = parse(Int, readline(io))
        readline(io)
        dim = parse(Int, readline(io))
        centers = zeros(dim, ncell)
        radii = zeros(1, ncell)
        readline(io)
        for icell = 1:ncell
            tmp = parse.(Float64, split(readline(io)))
            centers[:, icell] = tmp[2:1+dim]
            radii[icell] = tmp[2+dim]
        end
        (; centers, radii)
    end
end
