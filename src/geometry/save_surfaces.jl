"""
    save_surfaces(filename, surfaces)

Save surface triangulation.
The files may be passed to Tetgen with `filename`.node.
"""
function save_surfaces(filename, surfaces)

    default_refinement = 0.1

    # Extract surface triangulation
    (; points, facets, facetmarkers, regions) = surfaces

    npoint = size(points, 2)
    nfacet = size(facets, 2)
    nregion = size(regions, 2)

    # Write list of nodes in .node file
    open(filename * ".node", "w") do io
        write(io, "# Part 1 - node list\n")
        write(io, "# node count, 3D, no attribute, no boundary marker\n")
        write(io, "$npoint 3 0 0\n")
        write(io, "# Node index, node coordinates\n")
        for i = 1:npoint
            write(io, @sprintf "%d %26.16f %26.16f %26.16f\n" i points[:, i]...)
        end
    end

    # Write .poly file
    open(filename * ".poly", "w") do io

        # Write list of holes (refer to separate file)
        write(io, "# Part 1 - node list\n")
        write(io, "#  0 indicates the node list is stored in file .node\n")
        write(io, "0\n")

        # Write list of facets
        write(io, "# Part 2 - facet list\n")
        write(io, "# facet count, yes boundary marker\n")
        write(io, "$nfacet 1\n")
        write(io, "# Node index, node coordinates\n")
        for ifacet = 1:nfacet
            write(io, "1 0 $(facetmarkers[ifacet])\n")
            write(io, @sprintf "%d %d %d %d \n" 3 facets[:, ifacet]...)
        end

        # Write list of holes (empty)
        write(io, "# Part 3 - hole list\n")
        write(io, "0\n")

        # Write list of interior points with their corresponding compartment
        write(io, "# Part 4 - region list\n")
        write(io, "$nregion\n")
        for i = 1:nregion
            write(
                io,
                @sprintf "%d %f %f %f %d %f\n" i regions[:, i]... i default_refinement
            )
        end

    end

end
