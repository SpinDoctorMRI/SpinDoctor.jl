"""
    save_mesh(mesh_all, filename)

Save mesh in the Tetgen output format.
"""
function save_mesh(mesh_all, filename)
    # Extract mesh
    (; points, facets, elements, facetmarkers, elementmarkers) = mesh_all

    @info "Saving mesh in Tetgen output format at " * filename

    # Save points
    dim, npoint = size(points)
    open(filename * ".node", "w") do io
        write(io, "$npoint $dim 0 0\n")
        for ipoint = 1:npoint
            if dim == 2
                write(io, @sprintf("%d %f %f\n", ipoint, points[:, ipoint]...))
            elseif dim == 3
                write(io, @sprintf("%d %f %f %f\n", ipoint, points[:, ipoint]...))
            end
        end
    end

    # Read facets and their associated boundaries
    nfacet = size(facets, 2)
    open(filename * ".face", "w") do io
        write(io, "$nfacet 1\n")
        for i = 1:nfacet
            if dim == 2
                write(io, @sprintf("%d %d %d %d\n", i, facets[:, i]..., facetmarkers[i]))
            elseif dim == 3
                write(io, @sprintf("%d %d %d %d %d\n", i, facets[:, i]..., facetmarkers[i]))
            end
        end
    end

    # Read elements and their associated compartments
    nelement = size(elements, 2)
    open(filename * ".ele", "w") do io
        write(io, "$nelement $(dim+1) 1\n")
        for i = 1:nelement
            if dim == 2
                write(
                    io,
                    @sprintf("%d %d %d %d %d\n", i, elements[:, i]..., elementmarkers[i])
                )
            elseif dim == 3
                write(
                    io,
                    @sprintf(
                        "%d %d %d %d %d %d\n",
                        i,
                        elements[:, i]...,
                        elementmarkers[i]
                    )
                )
            end
        end
    end

end
