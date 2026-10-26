# Reading and writing geometries in text formats: Gmsh volume meshes, and the files of the
# MATLAB version of SpinDoctor (cell configurations, and surface and volume meshes in the
# TetGen format)

"""
    points, elements = read_gmsh(filename)

Read the points and tetrahedra of a volume mesh in the Gmsh format (version 2, as in the
RealNeuronMeshes repository). Other elements (e.g. surface triangles) are skipped.
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

        # Elements: number, type (4 for tetrahedra), number of tags, tags, and points
        readline(io)
        nelement = parse(Int, readline(io))
        elements = Int[]
        for _ = 1:nelement
            fields = parse.(Int, split(readline(io)))
            fields[2] == 4 && append!(elements, fields[(4+fields[3]):end])
        end

        points, reshape(elements, 4, :)
    end
end

"""
    read_cells(filename)

Read a cell configuration `(; centers, radii)`, with the fields `box` and `gap` if present
in the file.
"""
function read_cells(filename)
    @info "Loading cell configuration from " * filename
    open(filename, "r") do io
        readline(io)
        ncell = parse(Int, readline(io))
        readline(io)
        dim = parse(Int, readline(io))
        hasradii = occursin("Radius", readline(io))
        centers = zeros(dim, ncell)
        radii = zeros(1, ncell)
        for icell = 1:ncell
            tmp = parse.(Float64, split(readline(io)))
            centers[:, icell] = tmp[2:(1+dim)]
            hasradii && (radii[icell] = tmp[2+dim])
        end
        cells = hasradii ? (; centers, radii) : (; centers)
        while !eof(io)
            line = readline(io)
            if line == "Box (lower, upper):"
                lower = parse.(Float64, split(readline(io)))
                upper = parse.(Float64, split(readline(io)))
                cells = (; cells..., box = (lower, upper))
            elseif line == "Gap:"
                cells = (; cells..., gap = parse(Float64, readline(io)))
            end
        end
        cells
    end
end

"""
    save_cells(filename, cells)

Save a cell configuration `(; centers, radii)`, and its fields `box` and `gap` if present
(e.g. from [`periodic_packing`](@ref) or [`VoronoiSetup`](@ref), whose cells have no
radii). The file contains the spatial dimension: 2 for disks and cylinders, 3 for spheres.
"""
function save_cells(filename, cells)
    (; centers) = cells
    d, ncell = size(centers)
    hasradii = haskey(cells, :radii)
    @info "Writing cell geometry to " * filename
    open(filename, "w") do io
        println(io, "Number of cells:")
        println(io, ncell)
        println(io, "Dimension:")
        println(io, d)
        coordinates = d == 2 ? "x,y" : "x,y,z"
        println(io, "Number n, Center ($coordinates)", hasradii ? ", Radius r" : "")
        for i = 1:ncell
            values = hasradii ? [centers[:, i]; cells.radii[i]] : centers[:, i]
            println(io, i, " ", join((@sprintf("%.17g", v) for v in values), " "))
        end
        if haskey(cells, :gap)
            println(io, "Gap:")
            println(io, @sprintf("%.17g", cells.gap))
        end
        if haskey(cells, :box)
            println(io, "Box (lower, upper):")
            for corner in cells.box
                println(io, join((@sprintf("%.17g", v) for v in corner), " "))
            end
        end
    end
end

"""
    read_surfaces(filename)

Read a surface triangulation from the files `filename.node` and `filename.poly`.
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
    facets, facetmarkers, regions, regionmarkers = open(filename * ".poly", "r") do io

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
        regionmarkers = zeros(Int, nregion)
        for iregion = 1:nregion
            tmp = parse.(Float64, split(readline(io)))
            regions[:, iregion] = tmp[2:(2+dim-1)]
            regionmarkers[iregion] = round(Int, tmp[2+dim])
        end

        facets, facetmarkers, regions, regionmarkers
    end

    (; points, facets, facetmarkers, regions, regionmarkers)
end

"""
    save_surfaces(filename, surfaces)

Save a surface triangulation to the files `filename.node` and `filename.poly`, which may
be passed to TetGen.
"""
function save_surfaces(filename, surfaces)
    default_refinement = 0.1

    # Extract surface triangulation (polygon facets are split into triangles)
    (; points, regions) = surfaces
    facets, facetmarkers = triangulate_facets(surfaces)
    regionmarkers = get(surfaces, :regionmarkers, 1:size(regions, 2))

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
                        regionmarkers[i],
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
                        regionmarkers[i],
                        default_refinement
                    )
                )
            end
        end
    end
end

"""
    read_mesh(filename)

Read a global mesh from the TetGen output files `filename.node`, `filename.face` and
`filename.ele`.
"""
function read_mesh(filename)
    @info "Reading TetGen mesh from " * filename

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
    save_mesh(filename, mesh_all)

Save a global mesh in the TetGen output format (`filename.node`, `filename.face` and
`filename.ele`).
"""
function save_mesh(filename, mesh_all)
    # Extract mesh
    (; points, facets, elements, facetmarkers, elementmarkers) = mesh_all

    @info "Saving mesh in the TetGen format at " * filename

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

    # Save facets and their boundary numbers
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

    # Save elements and their compartment numbers
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
