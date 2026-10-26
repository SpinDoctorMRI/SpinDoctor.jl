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
            centers[:, icell] = tmp[2:(1+dim)]
            radii[icell] = tmp[2+dim]
        end
        (; centers, radii)
    end
end

"""
    save_cells(cells, cellfilename)

Write geometrical cell configuration to file.
The file contains the spatial dimension: 2 for cylinders, 3 for spheres.
"""
function save_cells(cells, cellfilename)

    # Extract cell parameters
    (; centers, radii) = cells

    # Sizes
    d, ncell = size(centers)

    # Save cell geometry to file
    @info "Writing cell geometry to " * cellfilename

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
                write(io, @sprintf "%d %.17g %.17g %.17g\n" i centers[:, i]... radii[i])
            else
                write(
                    io,
                    @sprintf "%d %.17g %.17g %.17g %.17g\n" i centers[:, i]... radii[i]
                )
            end
        end
    end
end

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
            points[:, ipoint] = tmp[2:(2+dim-1)]
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
            facets[:, ifacet] = tmp[2:(2+dim-1)]
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
            regions[:, iregion] = tmp[2:(2+dim-1)]
        end

        facets, facetmarkers, regions
    end

    (; points, facets, facetmarkers, regions)
end

"""
    save_surfaces(filename, surfaces)

Save surface triangulation.
The files may be passed to Tetgen with `filename`.node.
"""
function save_surfaces(filename, surfaces)
    default_refinement = 0.1

    # Extract surface triangulation
    (; points, facets, facetmarkers, regions) = surfaces

    dim, npoint = size(points)
    nfacet = size(facets, 2)
    nregion = size(regions, 2)

    # Write list of nodes in .node file
    open(filename * ".node", "w") do io
        write(io, "# Part 1 - node list\n")
        write(io, "# node count, dimension, no attribute, no boundary marker\n")
        write(io, "$npoint $dim 0 0\n")
        write(io, "# Node index, node coordinates\n")
        for i = 1:npoint
            if dim == 2
                write(io, @sprintf("%d %26.16f %26.16f\n", i, points[:, i]...))
            elseif dim == 3
                write(io, @sprintf("%d %26.16f %26.16f %26.16f\n", i, points[:, i]...))
            end
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
            if dim == 2
                write(io, @sprintf "%d %d %d \n" 2 facets[:, ifacet]...)
            elseif dim == 3
                write(io, @sprintf "%d %d %d %d \n" 3 facets[:, ifacet]...)
            end
        end

        # Write list of holes (empty)
        write(io, "# Part 3 - hole list\n")
        write(io, "0\n")

        # Write list of interior points with their corresponding compartment
        write(io, "# Part 4 - region list\n")
        write(io, "$nregion\n")
        for i = 1:nregion
            if dim == 2
                write(
                    io,
                    @sprintf(
                        "%d %.17g %.17g %d %f\n",
                        i,
                        regions[:, i]...,
                        i,
                        default_refinement
                    )
                )
            elseif dim == 3
                write(
                    io,
                    @sprintf(
                        "%d %.17g %.17g %.17g %d %f\n",
                        i,
                        regions[:, i]...,
                        i,
                        default_refinement
                    )
                )
            end
        end
    end
end

"""
    read_mesh(filename)
    
Read mesh from Tetgen format file.
"""
function read_mesh(filename)
    @info "Reading Tetgen FE mesh from " * filename

    # Read points
    points, dim = open(filename * ".node", "r") do io
        npoint, dim, = parse.(Int, split(readline(io)))
        points = zeros(dim, npoint)
        for ipoint = 1:npoint
            points[:, ipoint] = parse.(Float64, split(readline(io))[2:(1+dim)])
        end
        points, dim
    end

    # Read facets and their associated boundaries
    facets, facetmarkers = open(filename * ".face", "r") do io
        nfacet = parse(Int, split(readline(io))[1])
        facets = zeros(Int, dim, nfacet)
        facetmarkers = zeros(Int, nfacet)
        for ifacet = 1:nfacet
            tmp = parse.(Int, split(readline(io)))
            facets[:, ifacet] = tmp[2:(1+dim)]
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
            elements[:, ielement] = tmp[2:(2+dim)]
            elementmarkers[ielement] = tmp[3+dim]
        end
        elements, elementmarkers
    end

    (; points, facets, facetmarkers, elements, elementmarkers)
end

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
                write(io, @sprintf("%d %.17g %.17g\n", ipoint, points[:, ipoint]...))
            elseif dim == 3
                write(io, @sprintf("%d %.17g %.17g %.17g\n", ipoint, points[:, ipoint]...))
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
