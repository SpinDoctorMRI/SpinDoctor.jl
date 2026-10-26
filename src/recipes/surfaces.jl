"""
    create_surfaces(setup, cells)

Create surface triangulation from a setup and cells.

Returns a named tuple containing:

  - `points` (size `dim × npoint`),
  - `facets` (edges in 2D) (size `dim × nfacet`),
  - `facetmarkers` (edge markers in 2D) (size `nfacet`),
  - `regions` (one point inside each region) (size `dim × nregion`)
"""
function create_surfaces end

function create_surfaces(setup::PlateSetup{T}, _) where {T}
    (; widths, depth) = setup

    n = length(widths)
    w = sum(widths) / 2
    x = cumsum([-w; widths])
    y = depth / 2

    points = zeros(T, 2, 2(n + 1))
    for i = 1:(n+1)
        points[:, (1+2(i-1)):2i] = [
            x[i] x[i]
            -y y
        ]
    end

    edges = fill(0, 2, 1 + 3n)
    edges[:, 1] = [1, 2]
    for i = 1:n
        inds = 1 + 3(i - 1) .+ (1:3)
        edges[:, inds] .= 2(i - 1) .+ [
            1 4 3
            3 2 4
        ]
    end

    ninterface = (n - 1) * n ÷ 2
    edgemarkers = fill(0, 1 + 3n)

    # Left and first side walls
    edgemarkers[1:3] .= ninterface + 1

    # Interfaces and side walls after interface
    for i = 1:(n-1)
        edgemarkers[3i+1] = (i - 1) * (2 * n - i) ÷ 2 + 1
        edgemarkers[3i+2] = ninterface + 1 + i
        edgemarkers[3i+3] = ninterface + 1 + i
    end

    # Right
    edgemarkers[end] = ninterface + n

    regions = zeros(T, 2, n)
    regions[1, :] .= (x[1:(end-1)] .+ x[2:end]) ./ 2

    (; points, facets = edges, facetmarkers = edgemarkers, regions)
end

function create_surfaces(setup::Union{DiskSetup,SphereSetup}, cells)
    (; radii, centers) = cells
    (; ncell, layersizes, rmin, rmax, nsidewall, ecs) = setup

    dim = size(centers, 1)

    # Minimum number of side facets
    if dim == 2
        nside_min = 8
    elseif dim == 3
        nside_min = 20
    end

    nlayer = length(layersizes)
    rmean = (rmin + rmax) / 2

    # Create empty arrays
    cell_points = [zeros(0, 0) for _ = 1:ncell, _ = 1:nlayer]
    cell_facets = [zeros(Int, 0, 0) for _ = 1:ncell, _ = 1:nlayer]
    cell_regions = [zeros(0) for _ = 1:ncell, _ = 1:nlayer]
    for i = 1:nlayer, j = 1:ncell
        r = layersizes[i] * radii[j]
        if dim == 2
            n = max(nside_min, round(Int, nsidewall * r / rmean))
            angles = 2π * collect(0:(n-1)) / n
            cell_points[j, i] = r * [cos.(angles'); sin.(angles')] .+ centers[:, j]
            cell_facets[j, i] = [(1:(n-1))' n; (2:n)' 1]
        elseif dim == 3
            n = max(nside_min, round(Int, nsidewall * (r / rmean)^2))
            cell_points[j, i] = centers[:, j] .+ r .* create_fibonacci_sphere(n)
            cell_facets[j, i] = convexhull(cell_points[j, i])[1]
        end
        if i == 1
            cell_regions[j, i] = centers[:, j]
        else
            dir = zeros(dim)
            dir[1] = 1
            cell_regions[j, i] =
                centers[:, j] + dir * (layersizes[i-1] + layersizes[i]) / 2 * radii[j]
        end
    end

    # Cell points
    points = reduce(hcat, cell_points)

    # Create ECS surface
    ecs_surface = create_surfaces(ecs, points)
    points = [points ecs_surface.points]

    # Regions
    regions = [reduce(hcat, cell_regions) ecs_surface.region]
    nregion = size(regions, 2)

    # Take into account previous points in edges
    nprevious = 0
    for i = 1:nlayer, j = 1:ncell
        cell_facets[j, i] .+= nprevious
        nprevious += size(cell_points[j, i], 2)
    end
    ecs_surface.facets .+= nprevious
    facets = [reduce(hcat, reshape(cell_facets, :)) ecs_surface.facets]

    # The boundaries are ordered as follows:
    #   Number    Boundary
    #      1       Γ_1_2 (we skip Γ_2_1)
    #      2       Γ_1_3 (we skip Γ_3_1)
    #      ⋮         ⋮
    #  (n-1)n/2    Γ_n_n
    # (n-1)n/2+1    Γ_1
    #                ⋮
    # (n-1)n/2+n    Γ_n
    facetnumber(i, j) = (i - 1) * (2 * nregion - i) ÷ 2 + (j - i)
    facetmarkers = Int[]
    k = 1
    for i = 1:nlayer, j = 1:ncell
        if i < nlayer
            # (i, j) is linked to (i + 1, j)
            m = facetnumber(k, k + ncell)
        elseif !(ecs isa NoECS)
            # The outer layer is linked to the ECS (last region)
            m = facetnumber(k, nregion)
        else
            # The outer region k has an outer boundary Γ_k
            m = (nregion - 1) * nregion ÷ 2 + k
        end
        append!(facetmarkers, fill(m, size(cell_facets[j, i], 2)))
        k += 1
    end

    # ECS is last region, with outer boundary Γ_nregion
    append!(
        facetmarkers,
        fill((nregion - 1) * nregion ÷ 2 + nregion, size(ecs_surface.facets, 2)),
    )

    (; points, facets, facetmarkers, regions)
end

function create_surfaces(setup::ExtrusionSetup, cells)
    (; groundsetup, height, refinement) = setup

    # Create ground surface triangulation using 2D disk setup
    groundfacets = create_surfaces(groundsetup, cells)
    groundmesh = create_mesh(groundfacets, groundsetup.refinement)

    # groundfacets.facetmarkers
    # plot_surfaces(groundfacets, 1:10)
    # plot_mesh(split_mesh(groundmesh), 1:4)

    groundregions = groundfacets.regions
    nregion = size(groundregions, 2)

    # 2D facets are 3D edges, 2D elements are 3D facets
    points = groundmesh.points
    edges = groundmesh.facets
    edgemarkers = groundmesh.facetmarkers
    facets = groundmesh.elements
    facetmarkers = groundmesh.elementmarkers

    # Copy ground points to the two 3D planes top and bottom
    npoint = size(points, 2)
    z = fill(height / 2, 1, npoint)
    points = [points points; -z z]

    # Create two side facet triangles from each ground edge. They form a
    # rectangle representing the extruded edge. The top layer points are
    # numbered `npoint` higher than the bottom layer
    le = edges[1:1, :]
    ri = edges[2:2, :]
    sidefacets = [
        le le
        ri (ri .+ npoint)
        (ri .+ npoint) (le .+ npoint)
    ]

    # The boundaries are ordered as follows:
    #   Γ_1_2 (we skip Γ_2_1)
    #   Γ_1_3 (we skip Γ_3_1)
    #   ⋮
    #   Γ_nregion_nregion
    #   Γ_1
    #   ⋮
    #   Γ_nregion
    sidefacetmarkers = [edgemarkers; edgemarkers]

    # Bottom facets, top facets and side facets
    ninterface = (nregion - 1) * nregion ÷ 2
    facets = [sidefacets facets (facets .+ npoint)]
    facetmarkers =
        [sidefacetmarkers; facetmarkers .+ ninterface; facetmarkers .+ ninterface]

    # Put a point in the middle of each region (z = 0 is in the middle of
    # bottom and top layer)
    regions = [groundregions[1:2, :]; zeros(1, nregion)]

    (; points, facets, facetmarkers, regions)
end

function create_surfaces(setup::NeuronSetup, meshfile)
    (; ecs) = setup

    @info "Reading neuron FE mesh from " * meshfile
    points, elements = read_gmsh(meshfile)

    # Boundary facets: the element faces that belong to one element only
    count = Dict{SVector{3,Int},Int}()
    for e ∈ axes(elements, 2), skip = 1:4
        f = sort(SVector{3}(ntuple(k -> elements[k < skip ? k : k + 1, e], Val(3))))
        count[f] = get(count, f, 0) + 1
    end
    facets = reduce(hcat, [Vector(f) for (f, c) ∈ count if c == 1])

    volumes, centers = get_mesh_volumes(points, elements)
    centermass = centers * volumes / sum(volumes)

    facetmarkers = ones(Int, size(facets, 2))
    regions = centermass

    # Create ECS around points
    ecs_surface = create_surfaces(ecs, points)

    npoint = size(points, 2)
    points = [points ecs_surface.points]
    facets = [facets ecs_surface.facets .+ npoint]
    append!(facetmarkers, fill(2, size(ecs_surface.facets, 2)))
    regions = [regions ecs_surface.region]

    (; points, facets, facetmarkers, regions)
end

function create_surfaces(::NoECS, points)
    T = eltype(points)
    dim = size(points, 1)
    (; points = zeros(T, dim, 0), facets = zeros(Int, dim, 0), region = zeros(T, dim, 0))
end

function create_surfaces(ecs::BoxECS{T}, points) where {T}
    dim = size(points, 1)
    (; margin) = ecs

    # Determine bounds of domain
    pmin = minimum(points; dims = 2) .- margin
    pmax = maximum(points; dims = 2) .+ margin

    if dim == 2
        # Define four corners of box ECS and their corresponding edges
        points = [
            pmin[1] pmax[1] pmax[1] pmin[1]
            pmin[2] pmin[2] pmax[2] pmax[2]
        ]
        facets = [
            1 2 3 4
            2 3 4 1
        ]
    elseif dim == 3
        # Define eight corners of box ECS and their corresponding edges
        points = [
            pmin[1] pmax[1] pmax[1] pmin[1] pmin[1] pmax[1] pmax[1] pmin[1]
            pmin[2] pmin[2] pmax[2] pmax[2] pmin[2] pmin[2] pmax[2] pmax[2]
            pmin[3] pmin[3] pmin[3] pmin[3] pmax[3] pmax[3] pmax[3] pmax[3]
        ]
        facets = [
            1 1 1 1 2 2 3 3 4 4 5 5
            2 3 2 6 3 7 4 8 1 5 6 7
            3 4 6 5 7 6 8 7 5 8 7 8
        ]
    end

    # The point half a margin in from a corner is guaranteed to be inside ECS
    region = pmin .+ margin / 2

    (; points, facets, region)
end

function create_surfaces(ecs::ConvexHullECS, points)
    (; margin, threshold, maxiter) = ecs
    threshold_inv = 1 / threshold
    dim = size(points, 1)

    facets, points = convexhull(points)
    nfacet = size(facets, 2)

    # Average point, is inside the convex hull by definition (convex
    # combination of points)
    c = sum(points; dims = 2) / size(points, 2)

    # Get outward pointing normals
    _, facet_centers, normals = get_mesh_surface(points, facets)
    o = sum((c .- facet_centers) .* normals; dims = 1)
    orientations = 1 .- 2 .* (o .> 0)
    normals .*= orientations

    # Add a point in the middle of each degenerate facet, and move it slightly
    # outwards
    if dim == 3
        for iter = 1:maxiter
            # Check for degenerate facets
            iadd = fill(false, nfacet)
            for i = 1:nfacet
                f = facets[:, i]
                a = points[:, f[2]] - points[:, f[1]]
                b = points[:, f[3]] - points[:, f[1]]
                c = points[:, f[3]] - points[:, f[2]]
                α = acos(dot(a, b) / norm(a) / norm(b))
                β = acos(dot(a, c) / norm(a) / norm(c))
                if abs(α) < threshold || abs(β) < threshold
                    iadd[i] = true
                end
            end

            padd = facet_centers[:, iadd] .+ 0.05 .* margin .* normals[:, iadd]
            points = [points padd]

            facets, points = convexhull(points)
            nfacet = size(facets, 2)

            # Average point, is inside the convex hull by definition (convex
            # combination of points)
            c = sum(points; dims = 2) / size(points, 2)

            # Get outward pointing normals
            _, facet_centers, normals = get_mesh_surface(points, facets)
            o = sum((c .- facet_centers) .* normals; dims = 1)
            orientations = 1 .- 2 .* (o .> 0)
            normals .*= orientations

            @info "Breaking up degenerate ECS facets" iter sum(iadd)
            sum(iadd) > 0 || break
        end
    end

    # Move each facet one margin outwards
    dirs = zero(points)
    for i = 1:size(facets, 2)
        # dirs[:, facets[:, i]] .+= normals[:, i]
        @views dirs[:, facets[:, i]] .+= normals[:, i]
    end
    ihull = map(c -> any(!≈(0), c), eachcol(dirs))
    points[:, ihull] .+=
        margin .* dirs[:, ihull] ./ sqrt.(sum(abs2, dirs[:, ihull]; dims = 1))

    # The point half a margin out from the first point is guaranteed inside
    # the ECS
    i = findfirst(ihull)
    region = points[:, i] - margin / 2 * dirs[:, i] / norm(dirs[:, i])

    (; points, facets, region)
end
