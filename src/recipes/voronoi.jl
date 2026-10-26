# Polyhedral cells: shrunk Voronoi cells of seed points in a periodic box

"""
    VoronoiSetup(;
        ncell = 20,
        boxsize = 20.0,
        fraction = 0.8,
        lloyd = 10,
        meshsize = boxsize / ncell^(1 / 3) / 8,
        refinement = Inf,
    )

Setup recipe for `ncell` polyhedral cells in a periodic cubic box of side `boxsize`: the
Voronoi cells of random seed points (made more regular by `lloyd` Lloyd iterations),
separated by a uniform gap such that the cells fill the volume fraction `fraction`. The
cells are cut by the faces of the box, their periodic images enter on the opposite faces,
and the extracellular space fills the gaps (as for a [`PeriodicBoxECS`](@ref)).

Unlike packings of spheres, which cannot fill more than about 64 % of the volume,
polyhedral cells can fill any fraction below one. The surfaces are triangulated with
triangles of size `meshsize`, such that the opposite faces of the box have the same
triangles, and the tetrahedra have a volume of at most `refinement`. The mesh is periodic,
so that periodic boundary conditions may be used with `assemble_matrices(model; periodic =
true)`. The coefficients are given as for a [`SphereSetup`](@ref) with one layer.
"""
@with_kw struct VoronoiSetup{T} <: AbstractSetup{T}
    ncell::Int = 20
    @assert ncell > 0
    boxsize::T = 20.0
    @assert boxsize > 0
    fraction::T = 0.8
    @assert 0 < fraction < 1
    lloyd::Int = 10
    @assert lloyd >= 0
    meshsize::T = boxsize / ncell^(1 / 3) / 8
    @assert meshsize > 0
    refinement::T = Inf
    @assert refinement > 0
end

"""
    create_cells(setup::VoronoiSetup; rng = Random.default_rng())

Seed points (`centers`) of the Voronoi cells, after Lloyd iterations, the `gap` between
the cells giving the cell volume fraction, and the `box`. The seeds are translated such
that the box faces cut the cells away from their vertices and edges, so that the cut cells
have no slivers.
"""
function create_cells(setup::VoronoiSetup{T}; rng = Random.default_rng()) where {T}
    (; ncell, boxsize, fraction, lloyd) = setup
    L = T(boxsize)
    seeds = rand(rng, T, 3, ncell) .* L

    # Lloyd iterations: move the seeds to the centroids of their cells
    for _ = 1:lloyd
        cells = periodic_voronoi(seeds, L, zero(T))
        seeds = mod.(reduce(hcat, polyhedron_centroid.(cells)), L)
    end

    # Gap giving the volume fraction (the cell volumes decrease with the gap)
    volume(gap) = sum(polyhedron_volume, periodic_voronoi(seeds, L, gap)) / L^3
    a, b = zero(T), L / ncell^(1 / 3)
    for _ = 1:60
        gap = (a + b) / 2
        volume(gap) > fraction ? (a = gap) : (b = gap)
    end
    gap = (a + b) / 2

    # Translation of the cells maximizing the shortest edge created by the box faces
    lower, upper = zero(SVector{3,T}), fill(L, SVector{3,T})
    polyhedra = periodic_voronoi(seeds, L, gap)
    tol = sqrt(eps(T)) * L
    onbox(p) = any(@. abs(p - lower) <= tol || abs(p - upper) <= tol)
    function shortest_cut_edge(t)
        moved = [[[v + t for v in face] for face in polyhedron] for polyhedron in polyhedra]
        minimum(box_pieces(moved, lower, upper); init = L) do (_, piece)
            minimum(piece) do face
                minimum(eachindex(face)) do k
                    a, b = face[k], face[mod1(k + 1, length(face))]
                    onbox(a) || onbox(b) ? norm(b - a) : L
                end
            end
        end
    end
    translations = [[zero(SVector{3,T})]; [L * rand(rng, SVector{3,T}) for _ = 1:31]]
    t = argmax(shortest_cut_edge, translations)
    seeds = mod.(seeds .+ t, L)

    (; centers = seeds, gap, box = (zeros(T, 3), fill(L, 3)))
end

# Convex polyhedra are lists of faces, each a list of vertices, counterclockwise as seen
# from outside

# Voronoi cells of the seeds in a periodic cube of side L, with the faces between
# neighboring cells moved inwards by gap / 2. Each cell surrounds its seed.
function periodic_voronoi(seeds, L, gap)
    T = eltype(seeds)
    n = size(seeds, 2)
    shifts = [SVector{3,T}(s) * L for s in Iterators.product(-1:1, -1:1, -1:1)]
    map(1:n) do i
        c = SVector{3,T}(seeds[:, i])
        cell = box_polyhedron(c .- L / 2, c .+ L / 2)
        neighbors = [
            SVector{3,T}(seeds[:, j]) + s for j = 1:n, s in shifts if !(j == i && iszero(s))
        ]
        sort!(neighbors; by = q -> norm(q - c))
        for q in neighbors
            # Neighbors further than twice the cell radius do not cut the cell
            radius = maximum(v -> norm(v - c), Iterators.flatten(cell))
            norm(q - c) > 2 * radius + gap && break
            u = (q - c) / norm(q - c)
            cell = clip_polyhedron(cell, u, dot(u, (c + q) / 2) - gap / 2)
        end
        cell
    end
end

function box_polyhedron(lower, upper)
    x(i, j, k) =
        SVector((lower[1], upper[1])[i], (lower[2], upper[2])[j], (lower[3], upper[3])[k])
    [
        [x(1, 1, 1), x(1, 2, 1), x(2, 2, 1), x(2, 1, 1)], # z = lower
        [x(1, 1, 2), x(2, 1, 2), x(2, 2, 2), x(1, 2, 2)], # z = upper
        [x(1, 1, 1), x(2, 1, 1), x(2, 1, 2), x(1, 1, 2)], # y = lower
        [x(1, 2, 1), x(1, 2, 2), x(2, 2, 2), x(2, 2, 1)], # y = upper
        [x(1, 1, 1), x(1, 1, 2), x(1, 2, 2), x(1, 2, 1)], # x = lower
        [x(2, 1, 1), x(2, 2, 1), x(2, 2, 2), x(2, 1, 2)], # x = upper
    ]
end

# Keep the part of the convex polyhedron with u . x <= d
function clip_polyhedron(faces, u, d)
    isempty(faces) && return faces
    T = eltype(eltype(eltype(faces)))
    tol = sqrt(eps(T)) * maximum(v -> norm(v), Iterators.flatten(faces))
    clipped = Vector{Vector{SVector{3,T}}}()
    cap = SVector{3,T}[]
    for face in faces
        newface = SVector{3,T}[]
        m = length(face)
        for k = 1:m
            a, b = face[k], face[mod1(k + 1, m)]
            da, db = dot(u, a) - d, dot(u, b) - d
            da <= 0 && push!(newface, a)
            if (da < 0 && db > 0) || (da > 0 && db < 0)
                p = a + da / (da - db) * (b - a)
                push!(newface, p)
                push!(cap, p)
            elseif da == 0
                push!(cap, a)
            end
        end
        newface = dedup_polygon(newface, tol)
        length(newface) >= 3 && push!(clipped, newface)
    end

    # The new face in the plane, ordered counterclockwise around u
    cap = dedup_points(cap, tol)
    if length(cap) >= 3
        c = sum(cap) / length(cap)
        e1 = normalize(cap[1] - c)
        e2 = cross(u, e1)
        sort!(cap; by = p -> atan(dot(p - c, e2), dot(p - c, e1)))
        push!(clipped, cap)
    end
    clipped
end

function dedup_points(points, tol)
    unique = empty(points)
    for p in points
        any(q -> norm(p - q) <= tol, unique) || push!(unique, p)
    end
    unique
end

# Remove consecutive duplicate vertices of a polygon
function dedup_polygon(polygon, tol)
    out = empty(polygon)
    for p in polygon
        (isempty(out) || norm(p - out[end]) > tol) && push!(out, p)
    end
    length(out) > 1 && norm(out[1] - out[end]) <= tol && pop!(out)
    out
end

# The parts of the polyhedra and their periodic images inside the box, with the index of
# their polyhedron
function box_pieces(polyhedra, lower, upper)
    T = eltype(lower)
    L = upper - lower
    tol = sqrt(eps(T)) * norm(L)
    pieces = Tuple{Int,Vector{Vector{SVector{3,T}}}}[]
    for (i, polyhedron) in enumerate(polyhedra), s in Iterators.product(-1:1, -1:1, -1:1)
        piece = [[v + SVector{3,T}(s) .* L for v in face] for face in polyhedron]
        for d = 1:3
            piece = clip_polyhedron(piece, -SVector{3,T}(1:3 .== d), -lower[d])
            piece = clip_polyhedron(piece, SVector{3,T}(1:3 .== d), upper[d])
            length(piece) >= 4 || break
        end
        length(piece) >= 4 && polyhedron_volume(piece) > tol^2 && push!(pieces, (i, piece))
    end
    pieces
end

# Unit normal of a planar polygon (Newell's method), outward for counterclockwise vertices
# as seen from outside
function polygon_normal(polygon)
    n = sum(eachindex(polygon)) do k
        a, b = polygon[k], polygon[mod1(k + 1, length(polygon))]
        cross(a, b)
    end
    normalize(n)
end

# Whether the point p is inside the convex polyhedron
inpolyhedron(p, faces) =
    all(face -> dot(polygon_normal(face), p - sum(face) / length(face)) < 0, faces)

function polyhedron_volume(faces)
    c = sum(sum, faces) / sum(length, faces)
    sum(faces) do face
        sum(2:(length(face)-1)) do k
            abs(dot(face[1] - c, cross(face[k] - c, face[k+1] - c))) / 6
        end
    end
end

function polyhedron_centroid(faces)
    c = sum(sum, faces) / sum(length, faces)
    V = zero(eltype(c))
    m = zero(c)
    for face in faces, k = 2:(length(face)-1)
        v = abs(dot(face[1] - c, cross(face[k] - c, face[k+1] - c))) / 6
        V += v
        m += v * (c + face[1] + face[k] + face[k+1]) / 4
    end
    m / V
end

# The cells (with their periodic images) cut by the faces of the box, and the ECS between
# them, with a triangulation of all the surfaces with triangles of size `meshsize`. The
# edges of the polygons are split into segments of at most `meshsize`, and the polygons
# are triangulated in their planes. The lower faces of the box are triangulated, and the
# upper faces are their translates, so that the surface triangulation is periodic.
function create_surfaces(setup::VoronoiSetup{T}, cells) where {T}
    (; ncell, meshsize) = setup
    (; centers, gap, box) = cells
    lower, upper = SVector{3,T}(box[1]), SVector{3,T}(box[2])
    L = upper - lower
    nregion = ncell + 1
    tol = sqrt(eps(T)) * norm(L)
    maxarea = meshsize^2 / 2

    polyhedra = periodic_voronoi(centers, L[1], gap)
    pieces = box_pieces(polyhedra, lower, upper)

    points = SVector{3,T}[]
    function addpoint(p)
        # The vertices are shared between the faces of a piece, and with the box faces
        for (k, q) in enumerate(points)
            norm(p - q) <= tol && return k
        end
        push!(points, p)
        length(points)
    end

    # Face f of the box: axis div(f + 1, 2), lower for odd f
    faceaxis(f) = div(f + 1, 2)
    faceplane(f) = isodd(f) ? lower[faceaxis(f)] : upper[faceaxis(f)]
    onface(p, f) = abs(p[faceaxis(f)] - faceplane(f)) <= tol

    # Polygons of the cell pieces: membranes, and cell parts of the box faces
    regions = SVector{3,T}[]
    regionmarkers = Int[]
    membranes = Tuple{Int,Vector{Int}}[]
    boxpolygons = [Tuple{Int,Vector{Int}}[] for _ = 1:6]
    ecsseeds = SVector{3,T}[]
    for (i, piece) in pieces
        push!(regions, polyhedron_centroid(piece))
        push!(regionmarkers, i)
        for face in piece
            polygon = addpoint.(face)
            f = findfirst(f -> all(p -> onface(p, f), face), 1:6)
            if isnothing(f)
                push!(membranes, (i, polygon))
                # ECS point: halfway to the neighboring cell, in the gap
                p = sum(face) / length(face) + gap / 2 * polygon_normal(face)
                all(lower .< p .< upper) && push!(ecsseeds, p)
            else
                push!(boxpolygons[f], (i, polygon))
            end
        end
    end

    # Outlines of the box faces: the corners and the vertices on the edges of the box, in
    # counterclockwise order (as seen from the axis direction)
    corners = [lower + SVector{3,T}(i, j, k) .* L for i in (0, 1), j in (0, 1), k in (0, 1)]
    addpoint.(corners)
    outlines = map(1:6) do f
        d = faceaxis(f)
        e1, e2 = filter(!=(d), 1:3)
        onedge(p) =
            onface(p, f) && (
                abs(p[e1] - lower[e1]) <= tol ||
                abs(p[e1] - upper[e1]) <= tol ||
                abs(p[e2] - lower[e2]) <= tol ||
                abs(p[e2] - upper[e2]) <= tol
            )
        center = (lower + upper) / 2
        angle(k) = atan(points[k][e2] - center[e2], points[k][e1] - center[e1])
        sort(filter(k -> onedge(points[k]), eachindex(points)); by = angle)
    end

    # Split the edges into segments of at most `meshsize` (with the same points for the
    # edges of neighboring polygons)
    edgepoints = Dict{Tuple{Int,Int},Vector{Int}}()
    function chain(a, b)
        key = minmax(a, b)
        inner = get!(edgepoints, key) do
            p, q = points[key[1]], points[key[2]]
            nsplit = max(1, ceil(Int, norm(q - p) / meshsize))
            map(1:(nsplit-1)) do k
                push!(points, p + k / nsplit * (q - p))
                length(points)
            end
        end
        a < b ? inner : reverse(inner)
    end
    split(polygon) = reduce(
        vcat,
        [
            [polygon[k]; chain(polygon[k], polygon[mod1(k + 1, length(polygon))])] for
            k in eachindex(polygon)
        ],
    )
    outlines = split.(outlines)
    membranes = [(i, split(polygon)) for (i, polygon) in membranes]
    boxpolygons =
        [[(i, split(polygon)) for (i, polygon) in polygons] for polygons in boxpolygons]

    facets = Vector{Int}[]
    facetmarkers = Int[]

    # Triangulate a polygon (or a face with polygons inside) in the plane with the origin o
    # and the axes e1, e2, keeping the boundary points. Returns the triangles in global
    # numbering, the region attributes and the new points.
    function triangulate_plane(loops, seeds, attributes, o, e1, e2)
        local2global = unique(reduce(vcat, loops))
        global2local = Dict(k => i for (i, k) in enumerate(local2global))
        segments = unique([
            collect(
                minmax(
                    global2local[loop[k]],
                    global2local[loop[mod1(k + 1, length(loop))]],
                ),
            ) for loop in loops for k in eachindex(loop)
        ])
        coordinates(p) = [dot(p - o, e1), dot(p - o, e2)]
        triin = TriangulateIO()
        triin.pointlist = reduce(hcat, [coordinates(points[k]) for k in local2global])
        triin.segmentlist = reduce(hcat, segments)
        if !isempty(seeds)
            triin.regionlist = [
                reduce(hcat, coordinates.(seeds))
                attributes'
                fill(maxarea, 1, length(seeds))
            ]
        end
        triout = first(triangulate("pqYYa$(maxarea)AQ", triin))
        newpoints = [
            o + q[1] * e1 + q[2] * e2 for
            q in eachcol(triout.pointlist[:, (length(local2global)+1):end])
        ]
        triangles = map(eachcol(triout.trianglelist)) do t
            [
                k <= length(local2global) ? local2global[k] : -(k - length(local2global)) for k in t
            ]
        end
        attrs =
            isempty(seeds) ? zeros(Int, length(triangles)) :
            round.(Int, triout.triangleattributelist[:])
        triangles, attrs, newpoints
    end

    # Membranes
    for (i, polygon) in membranes
        vertices = points[polygon]
        n = polygon_normal(vertices)
        e1 = normalize(vertices[2] - vertices[1])
        triangles, _, newpoints =
            triangulate_plane([polygon], [], Int[], vertices[1], e1, cross(n, e1))
        offset = length(points)
        append!(points, newpoints)
        for t in triangles
            push!(facets, [k > 0 ? k : offset - k for k in t])
            push!(facetmarkers, interface_index(i, nregion, nregion))
        end
    end

    # Lower faces of the box, and their translates on the upper faces
    for f in (1, 3, 5)
        d = faceaxis(f)
        e1, e2 = SVector{3,T}.(eachcol(Matrix{T}(I, 3, 3)[:, filter(!=(d), 1:3)]))
        polygons = boxpolygons[f]
        seeds = [sum(points[polygon]) / length(polygon) for (_, polygon) in polygons]
        triangles, attrs, newpoints = triangulate_plane(
            [[outlines[f]]; last.(polygons)],
            seeds,
            first.(polygons),
            lower,
            e1,
            e2,
        )
        offset = length(points)
        append!(points, newpoints)
        triangles = [[k > 0 ? k : offset - k for k in t] for t in triangles]

        # The points of the upper face: the translates of the points on the lower face
        shift = L[d] * SVector{3,T}(1:3 .== d)
        lowerpoints = unique(reduce(vcat, triangles))
        upperpoints = findall(p -> onface(p, f + 1), points)
        matches = match_points(
            reduce(hcat, points[upperpoints]),
            reduce(hcat, points[lowerpoints] .+ Ref(shift)),
            tol,
        )
        translate = Dict{Int,Int}()
        for (k, m) in zip(lowerpoints, matches)
            if m > 0
                translate[k] = upperpoints[m]
            elseif k > offset
                push!(points, points[k] + shift)
                translate[k] = length(points)
            else
                error("The opposite faces of the box do not match")
            end
        end

        for (t, attr) in zip(triangles, attrs)
            marker = outer_boundary_index(attr == 0 ? nregion : attr, nregion)
            push!(facets, t, [translate[k] for k in t])
            push!(facetmarkers, marker, marker)
        end

        # ECS points slightly inside the box from the largest triangle of each connected
        # part of the ECS on the face. The ECS in the box may have several parts (e.g.
        # pockets at the edges of the box), which all touch the box faces.
        ecs = findall(iszero, attrs)
        parent = collect(eachindex(ecs))
        root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
        edgetriangle = Dict{Tuple{Int,Int},Int}()
        for (i, k) in enumerate(ecs), j = 1:3
            t = triangles[k]
            edge = minmax(t[j], t[mod1(j + 1, 3)])
            haskey(edgetriangle, edge) ? (parent[root(i)] = root(edgetriangle[edge])) :
            (edgetriangle[edge] = i)
        end
        inradius(t) = (
            (a, b, c) = points[t];
            norm(cross(b - a, c - a)) / (norm(b - a) + norm(c - b) + norm(a - c))
        )
        inward = SVector{3,T}(1:3 .== d)
        for r in unique(root.(eachindex(ecs)))
            k = ecs[argmax(
                i -> root(i) == r ? inradius(triangles[ecs[i]]) : -one(T),
                eachindex(ecs),
            )]
            p = sum(points[triangles[k]]) / 3
            dp = inradius(triangles[k]) / 1000 * inward
            push!(ecsseeds, p + dp, p + shift - dp)
        end
    end

    # Keep the ECS points outside the cells
    filter!(p -> !any(((_, piece),) -> inpolyhedron(p, piece), pieces), ecsseeds)
    isempty(ecsseeds) && error("No ECS found between the cells")
    append!(regions, ecsseeds)
    append!(regionmarkers, fill(nregion, length(ecsseeds)))

    (;
        points = reduce(hcat, points),
        facets = reduce(hcat, facets),
        facetmarkers,
        regions = reduce(hcat, regions),
        regionmarkers,
        box = (Vector(lower), Vector(upper)),
        meshsize,
    )
end

coefficients(setup::VoronoiSetup; kwargs...) =
    coefficients(SphereSetup(; ncell = setup.ncell, ecs = BoxECS()); kwargs...)
