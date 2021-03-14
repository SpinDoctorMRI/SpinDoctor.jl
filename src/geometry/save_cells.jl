"""
    save_cells(cells, cellfilename)

Write geometrical cell configuration to file.
The file contains the spatial dimension: 2 for cylinders, 3 for spheres.
"""
function save_cells(cells, cellfilename)

    # Extract cell parameters
    @unpack centers, radii = cells

    # Sizes
    d, ncell = size(centers)

    # Save cell geometry to file
    println("Writing cell geometry to " * cellfilename)

    open(cellfilename, "w") do io
        write(io, "Number of cells:\n")
        write(io, "$ncell\n")
        write(io, "Dimension:\n")
        write(io, "$d\n")
        c = "x,y"
        if d == 3
            c = c * ",z"
        end
        write(io, "Number n, Center ($c), Radius r\n")
        for i = 1:ncell
            if d == 2
                write(io, @sprintf "%d %g %g %g\n" i centers[:, i]... radii[i])
            else
                write(io, @sprintf "%d %g %g %g %g\n" i centers[:, i]... radii[i])
            end
        end
    end

end
