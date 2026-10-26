# Surface triangulations of the setup recipes and their extracellular spaces

"""
    create_surfaces(setup, cells)

Create the surface triangulation of a setup with the cell configuration `cells` (for
neurons: the file name of the volume mesh), to be meshed by [`create_mesh`](@ref).

Returns a named tuple containing:

- `points`: the points (size `dim`, `npoint`);
- `facets`: the triangles (edges in 2D) (size `dim`, `nfacet`);
- `facetmarkers`: the boundary number of each facet;
- `regions`: points inside the regions (size `dim`, `npoint`), by default one per region;
- `regionmarkers` (optional): the region number of each of these points, for regions
  with several parts;
- `box` and `meshsize` (samples of periodic tissue): the box `(lower, upper)`, whose
  opposite faces have the same points, and the size of the facets. The mesh keeps the
  facets on the box (in 2D) or all the facets (in 3D), so that the mesh is periodic.

The boundaries between `n` regions are numbered with the interfaces between regions
`i < j` first (row by row: 12, 13, ..., 1n, 23, ...), followed by the outer boundaries of
the regions `1, ..., n`.
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

    edgemarkers = fill(0, 1 + 3n)

    # Left and first side walls
    edgemarkers[1:3] .= outer_boundary_index(1, n)

    # Interfaces and side walls after interface
    for i = 1:(n-1)
        edgemarkers[3i+1] = interface_index(i, i + 1, n)
        edgemarkers[3i+2] = outer_boundary_index(i + 1, n)
        edgemarkers[3i+3] = outer_boundary_index(i + 1, n)
    end

    # Right
    edgemarkers[end] = outer_boundary_index(n, n)

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
            angles = 2pi * collect(0:(n-1)) / n
            cell_points[j, i] = r * [cos.(angles'); sin.(angles')] .+ centers[:, j]
            cell_facets[j, i] = [(1:(n-1))' n; (2:n)' 1]
        elseif dim == 3
            n = max(nside_min, round(Int, nsidewall * (r / rmean)^2))
            cell_points[j, i] = centers[:, j] .+ r .* fibonacci_sphere(n)
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

    # Boundary numbers: see `interface_index` and `outer_boundary_index`
    facetmarkers = Int[]
    k = 1
    for i = 1:nlayer, j = 1:ncell
        if i < nlayer
            # (i, j) is linked to (i + 1, j)
            m = interface_index(k, k + ncell, nregion)
        elseif !(ecs isa NoECS)
            # The outer layer is linked to the ECS (last region)
            m = interface_index(k, nregion, nregion)
        else
            # Outer boundary of the outer layer
            m = outer_boundary_index(k, nregion)
        end
        append!(facetmarkers, fill(m, size(cell_facets[j, i], 2)))
        k += 1
    end

    # The ECS is the last region
    append!(
        facetmarkers,
        fill(outer_boundary_index(nregion, nregion), size(ecs_surface.facets, 2)),
    )

    (; points, facets, facetmarkers, regions)
end

function create_surfaces(setup::ExtrusionSetup, cells)
    (; groundsetup, height, refinement) = setup

    # Create ground surface triangulation using 2D disk setup
    groundfacets = create_surfaces(groundsetup, cells)
    groundmesh = create_mesh(groundfacets, groundsetup.refinement)

    groundregions = groundfacets.regions
    regionmarkers = get(groundfacets, :regionmarkers, 1:size(groundregions, 2))
    nregion = maximum(regionmarkers)

    # 2D facets are 3D edges, 2D elements are 3D facets
    points = groundmesh.points
    edges = groundmesh.facets
    edgemarkers = groundmesh.facetmarkers
    facets = groundmesh.elements
    facetmarkers = groundmesh.elementmarkers

    # Copy the ground points to the layers between the bottom and the top. Periodic samples
    # are cut into layers of the mesh size, since the mesh keeps their facets, and the other
    # setups have one layer.
    periodic = haskey(groundfacets, :box)
    nlayer = periodic ? max(1, ceil(Int, height / groundfacets.meshsize)) : 1
    npoint = size(points, 2)
    z = range(-height / 2, height / 2, nlayer + 1)
    points = [repeat(points, 1, nlayer + 1); repeat(z', npoint)[:]']

    # Each ground edge is extruded to a rectangle in each layer, split into two triangles.
    # In periodic samples, the diagonal goes from the first point of the edge (in
    # lexicographic order of the coordinates) in the lower layer, so that the opposite faces
    # of the box have the same triangles.
    sidefacets = zeros(Int, 3, 0)
    if periodic
        for (k, (a, b)) in enumerate(eachcol(edges))
            if groundmesh.points[:, b] < groundmesh.points[:, a]
                edges[:, k] = [b, a]
            end
        end
    end
    for l = 1:nlayer
        a = edges[1:1, :] .+ (l - 1) * npoint
        b = edges[2:2, :] .+ (l - 1) * npoint
        sidefacets = [sidefacets [a a; b (b .+ npoint); (b .+ npoint) (a .+ npoint)]]
    end

    # Boundary numbers: see `interface_index` and `outer_boundary_index`
    sidefacetmarkers = repeat(edgemarkers, 2nlayer)

    # Side facets, bottom facets and top facets. The bottom and top facets of a region are
    # on its outer boundary.
    facets = [sidefacets facets (facets .+ nlayer * npoint)]
    top = outer_boundary_index.(facetmarkers, nregion)
    facetmarkers = [sidefacetmarkers; top; top]

    # Put the region points in the middle of the bottom and top layers (z = 0), or in the
    # middle of a layer
    zregion = isodd(nlayer) ? 0.0 : height / 2nlayer
    regions = [groundregions[1:2, :]; fill(zregion, 1, size(groundregions, 2))]

    surfaces = (; points, facets, facetmarkers, regions, regionmarkers)
    if periodic
        lower, upper = groundfacets.box
        box = ([lower; -height / 2], [upper; height / 2])
        surfaces = (; surfaces..., box, groundfacets.meshsize)
    end
    surfaces
end

function create_surfaces(setup::NeuronSetup, meshfile)
    (; ecs) = setup

    @info "Reading neuron FE mesh from " * meshfile
    points, elements = read_gmsh(meshfile)

    # Boundary facets: the element faces that belong to one element only
    count = Dict{SVector{3,Int},Int}()
    for e in axes(elements, 2), skip = 1:4
        f = sort(SVector{3}(ntuple(k -> elements[k < skip ? k : k + 1, e], Val(3))))
        count[f] = get(count, f, 0) + 1
    end
    facets = reduce(hcat, [Vector(f) for (f, c) in count if c == 1])

    volumes, centers = element_geometry(points, elements)
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
    _, facet_centers, normals = facet_geometry(points, facets)
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
                alpha = acos(dot(a, b) / norm(a) / norm(b))
                beta = acos(dot(a, c) / norm(a) / norm(c))
                if abs(alpha) < threshold || abs(beta) < threshold
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
            _, facet_centers, normals = facet_geometry(points, facets)
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
    ihull = map(!iszero, eachcol(dirs))
    points[:, ihull] .+=
        margin .* dirs[:, ihull] ./ sqrt.(sum(abs2, dirs[:, ihull]; dims = 1))

    # The point half a margin out from the first point is guaranteed inside
    # the ECS
    i = findfirst(ihull)
    region = points[:, i] - margin / 2 * dirs[:, i] / norm(dirs[:, i])

    (; points, facets, region)
end

# Disks in a periodic box: the disks (and their periodic images) are cut by the faces of
# the box, and the ECS fills the rest of the box
function create_surfaces(setup::DiskSetup{T,<:PeriodicBoxECS}, cells) where {T}
    (; ncell, layersizes, rmin, rmax, nsidewall) = setup
    (; centers, radii, box) = cells
    lower, upper = SVector{2,T}(box[1]), SVector{2,T}(box[2])
    L = upper - lower
    nlayer = length(layersizes)
    nregion = nlayer * ncell + 1
    rmean = (rmin + rmax) / 2
    tol = 1e-10 * norm(L)

    points = SVector{2,T}[]
    edges = SVector{2,Int}[]
    edgemarkers = Int[]
    regions = SVector{2,T}[]
    regionmarkers = Int[]
    boundarypoints = Int[] # Indices of the points on the box boundary
    pieces = Tuple{Int,Vector{SVector{2,T}}}[] # Region and polygon of each cell piece

    # Cut the polygons of all layers of all periodic images of the cells that reach into
    # the box
    for j = 1:ncell, shift in Iterators.product(-1:1, -1:1)
        c = SVector{2,T}(centers[:, j]) + SVector{2,T}(shift) .* L
        rout = layersizes[end] * radii[j]
        all(lower .- rout .< c .< upper .+ rout) || continue
        inner = SVector{2,T}[]
        for i = 1:nlayer
            r = layersizes[i] * radii[j]
            n = max(8, round(Int, nsidewall * r / rmean))
            # Vertices close to a box face are moved away from it (towards the center), and
            # cuts close to a box edge are moved onto it, so that the cuts leave no slivers
            # and no membrane lies on a face. An eighth of the distance between the
            # vertices and the chords of their neighbors keeps the polygons convex.
            snap = r * (1 - cos(2pi / n)) / 8
            polygon = [
                away_from_faces(
                    c + r * SVector(cos(2pi * k / n), sin(2pi * k / n)),
                    c,
                    box,
                    snap,
                ) for k = 0:(n-1)
            ]
            polygon = clip_polygon(polygon, lower, upper, snap, tol)
            if length(polygon) < 3
                inner = SVector{2,T}[]
                continue
            end
            region = (i - 1) * ncell + j
            push!(pieces, (region, polygon))

            # Edges: membranes inside the box (the edges on the box are added below)
            offset = length(points)
            append!(points, polygon)
            m = length(polygon)
            marker =
                i < nlayer ? interface_index(region, region + ncell, nregion) :
                interface_index(region, nregion, nregion)
            for k = 1:m
                a, b = polygon[k], polygon[mod1(k + 1, m)]
                onboundary(a, lower, upper, tol) && push!(boundarypoints, offset + k)
                if !onsameside(a, b, lower, upper, tol)
                    push!(edges, SVector(offset + k, offset + mod1(k + 1, m)))
                    push!(edgemarkers, marker)
                end
            end

            # Region points: inside the piece, and outside the layer inside
            if isempty(inner)
                # Convex piece
                push!(regions, sum(polygon) / m)
                push!(regionmarkers, region)
            else
                # Ring piece (possibly in several parts): move inwards from the middle of
                # each membrane edge
                rin = layersizes[i-1] * radii[j]
                for k = 1:m
                    a, b = polygon[k], polygon[mod1(k + 1, m)]
                    onsameside(a, b, lower, upper, tol) && continue
                    e = (a + b) / 2
                    for step in (r - rin) .* (1 / 2, 1 / 4, 1 / 8, 1 / 16)
                        q = e + step * (c - e) / norm(c - e)
                        if inpolygon(q, polygon) && !inpolygon(q, inner)
                            push!(regions, q)
                            push!(regionmarkers, region)
                            break
                        end
                    end
                end
            end
            inner = polygon
        end
    end

    # The box boundary, split at the points of the cut polygons
    corners = [lower, SVector(upper[1], lower[2]), upper, SVector(lower[1], upper[2])]
    append!(points, corners)
    append!(boundarypoints, (length(points)-3):length(points))
    perimeter(p) =
        if abs(p[2] - lower[2]) <= tol
            p[1] - lower[1]
        elseif abs(p[1] - upper[1]) <= tol
            L[1] + p[2] - lower[2]
        elseif abs(p[2] - upper[2]) <= tol
            L[1] + L[2] + upper[1] - p[1]
        else
            2L[1] + L[2] + upper[2] - p[2]
        end
    s = perimeter.(points[boundarypoints])
    p = sortperm(s)
    boundarypoints, s = boundarypoints[p], s[p]

    # Merge coinciding points on the boundary (e.g. polygon vertices at the corners)
    duplicate = [false; diff(s) .<= tol]
    rep = collect(1:length(points))
    for k in findall(duplicate)
        rep[boundarypoints[k]] = rep[boundarypoints[k-1]]
    end
    boundarypoints = boundarypoints[.!duplicate]
    edges = [SVector(rep[e[1]], rep[e[2]]) for e in edges]

    # The segments of the box boundary are split into pieces of at most the side length of
    # the polygon of mean radius. The opposite sides of the box have the same points, and
    # the mesh keeps them, so that the mesh is periodic.
    meshsize = 2pi * rmean / nsidewall
    nb = length(boundarypoints)
    for k = 1:nb
        a, b = boundarypoints[k], boundarypoints[mod1(k + 1, nb)]

        # The region inside the box next to the middle of the segment (the box is on the
        # left of the counterclockwise boundary)
        t = points[b] - points[a]
        m = (points[a] + points[b]) / 2 + SVector(-t[2], t[1]) / 1000
        inside = [region for (region, polygon) in pieces if inpolygon(m, polygon)]
        region = isempty(inside) ? nregion : minimum(inside)
        nsplit = max(1, ceil(Int, norm(t) / meshsize))
        chain = [a]
        for i = 1:(nsplit-1)
            push!(points, points[a] + i / nsplit * t)
            push!(rep, length(points))
            push!(chain, length(points))
        end
        push!(chain, b)
        for i = 1:nsplit
            push!(edges, SVector(chain[i], chain[i+1]))
            push!(edgemarkers, outer_boundary_index(region, nregion))
        end

        # Region point: the ECS and the layers may have several parts in the box (e.g.
        # between two cells and a face, or between an inner layer and a corner), which the
        # points next to the membranes do not all reach
        push!(regions, m)
        push!(regionmarkers, region)
    end
    nregion in regionmarkers || error("No ECS found on the boundary of the box")

    # Remove the merged points
    used = rep .== 1:length(points)
    newindex = cumsum(used)
    edges = [SVector(newindex[e[1]], newindex[e[2]]) for e in edges]
    points = points[used]

    (;
        points = reduce(hcat, points),
        facets = reduce(hcat, edges),
        facetmarkers = edgemarkers,
        regions = reduce(hcat, regions),
        regionmarkers,
        box = (Vector(lower), Vector(upper)),
        meshsize,
    )
end

create_surfaces(::PeriodicBoxECS, points) =
    error("PeriodicBoxECS is only available for disks and cylinders")

# Move the coordinates of the point p that are within `snap` of a face of the periodic box
# (or of one of its periodic images) to the distance `snap` from the face, on the side of
# the center c
away_from_faces(p, c, (lower, upper), snap) =
    map(p, c, SVector{length(p)}(lower), SVector{length(p)}(upper)) do x, xc, l, u
        y = l + (u - l) * round((x - l) / (u - l))
        abs(x - y) < snap ? y + (xc < y ? -snap : snap) : x
    end

# Move the coordinates of the point p that are within `snap` of a face of the periodic box
# (or of one of its periodic images) onto the face
snap_to_faces(p, (lower, upper), snap) =
    map(p, SVector{length(p)}(lower), SVector{length(p)}(upper)) do x, l, u
        k = round((x - l) / (u - l))
        y = k == 0 ? l : k == 1 ? u : l + k * (u - l)
        abs(x - y) < snap ? y : x
    end

# Clip a convex polygon to the box `[lower, upper]` (Sutherland-Hodgman), with the cuts
# within `snap` of the other faces moved onto them, removing points closer than `tol` to
# their predecessor
function clip_polygon(polygon, lower, upper, snap, tol)
    for d = 1:2, (bound, sign) in ((lower[d], 1), (upper[d], -1))
        isempty(polygon) && break
        inside(p) = sign * (p[d] - bound) >= 0
        clipped = empty(polygon)
        m = length(polygon)
        for k = 1:m
            a, b = polygon[k], polygon[mod1(k + 1, m)]
            inside(a) && push!(clipped, a)
            if inside(a) != inside(b)
                t = (bound - a[d]) / (b[d] - a[d])
                q = a + t * (b - a)
                q = snap_to_faces(Base.setindex(q, bound, d), (lower, upper), snap)
                push!(clipped, q)
            end
        end
        polygon = empty(clipped)
        for p in clipped
            (isempty(polygon) || norm(p - polygon[end]) > tol) && push!(polygon, p)
        end
        length(polygon) > 1 && norm(polygon[1] - polygon[end]) <= tol && pop!(polygon)
    end
    polygon
end

onboundary(p, lower, upper, tol) = any(@. abs(p - lower) <= tol || abs(p - upper) <= tol)

# Whether the points a and b lie on the same face of the box
onsameside(a, b, lower, upper, tol) = any(
    d ->
        (abs(a[d] - lower[d]) <= tol && abs(b[d] - lower[d]) <= tol) ||
        (abs(a[d] - upper[d]) <= tol && abs(b[d] - upper[d]) <= tol),
    1:length(a),
)

# Whether the point p is inside the counterclockwise convex polygon
function inpolygon(p, polygon)
    m = length(polygon)
    all(1:m) do k
        a, b = polygon[k], polygon[mod1(k + 1, m)]
        (b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1]) > 0
    end
end

"""
    facets, facetmarkers = triangulate_facets(surfaces)

The facets of a surface triangulation as a matrix of triangles (edges in 2D), with polygon
facets (e.g. of [`VoronoiSetup`](@ref)) split into triangles.
"""
function triangulate_facets(surfaces)
    (; facets, facetmarkers) = surfaces
    facets isa AbstractMatrix && return facets, facetmarkers
    triangles = Vector{Int}[]
    markers = Int[]
    for (f, marker) in zip(facets, facetmarkers)
        for k = 2:(length(f)-1)
            push!(triangles, [f[1], f[k], f[k+1]])
            push!(markers, marker)
        end
    end
    reduce(hcat, triangles), markers
end
