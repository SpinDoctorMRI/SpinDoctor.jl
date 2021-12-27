"""
    save_tetgen(mesh_all, filename)
    
Save mesh in the Tetgen output format.
"""
function save_tetgen(mesh_all, filename)

    # Extract mesh
    @unpack points, facets, elements, facetmarkers, elementmarkers = mesh_all

    @info "Saving mesh in Tetgen output format at " * filename

    # Save points
    npoint = size(points, 2)
    open(filename * ".node", "w") do io
        write(io, "$npoint 3 0 0\n")
        for ipoint = 1:npoint
            write(io, @sprintf "%d %f %f %f\n" ipoint points[:, ipoint]...)
        end
    end

    # Read facets and their associated boundaries
    nfacet = size(facets, 2)
    open(filename * ".face", "w") do io
        write(io, "$nfacet 1\n")
        for i = 1:nfacet
            write(io, @sprintf "%d %d %d %d %d\n" i facets[:, i]... facetmarkers[i])
        end
    end

    # Read elements and their associated compartments
    nelement = size(elements, 2)
    open(filename * ".ele", "w") do io
        write(io, "$nelement 4 1\n")
        for i = 1:nelement
            write(io, @sprintf "%d %d %d %d %d %d\n" i elements[:, i]... elementmarkers[i])
        end
    end

end
