# Polyhedral cells: shrunk Voronoi cells of seed points in a periodic box

"""
    VoronoiSetup(;
        ncell = 20,
        boxsize = 20.0,
        fraction = 0.8,
        lloyd = 10,
        refinement = Inf,
    )

Setup recipe for `ncell` polyhedral cells in a periodic cubic box of side `boxsize`: the
Voronoi cells of random seed points (made more regular by `lloyd` Lloyd iterations),
separated by a uniform gap such that the cells fill the volume fraction `fraction`. The
cells are cut by the faces of the box, their periodic images enter on the opposite faces,
and the extracellular space fills the gaps (as for a [`PeriodicBoxECS`](@ref)).

Unlike packings of spheres, which cannot fill more than about 64 % of the volume,
polyhedral cells can fill any fraction below one. The tetrahedra have a volume of at most
`refinement`. The coefficients are given as for a [`SphereSetup`](@ref) with one layer.
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
# them. The facets are polygons, and the parts of the box faces in the ECS are triangulated.
function create_surfaces(setup::VoronoiSetup{T}, cells) where {T}
    (; ncell) = setup
    (; centers, gap, box) = cells
    lower, upper = SVector{3,T}(box[1]), SVector{3,T}(box[2])
    L = upper - lower
    nregion = ncell + 1
    tol = sqrt(eps(T)) * norm(L)

    polyhedra = periodic_voronoi(centers, L[1], gap)

    points = SVector{3,T}[]
    function addpoint(p)
        # Points are shared between the faces of a piece, and with the box faces
        for (k, q) in enumerate(points)
            norm(p - q) <= tol && return k
        end
        push!(points, p)
        length(points)
    end

    facets = Any[]
    facetmarkers = Int[]
    regions = SVector{3,T}[]
    regionmarkers = Int[]
    boxpolygons = [Vector{Int}[] for _ = 1:6] # Cell polygons on each box face
    onface(p, f) = abs(p[div(f + 1, 2)] - (isodd(f) ? lower : upper)[div(f + 1, 2)]) <= tol
    ecsseeds = SVector{3,T}[]

    pieces = box_pieces(polyhedra, lower, upper)
    for (i, piece) in pieces
        push!(regions, polyhedron_centroid(piece))
        push!(regionmarkers, i)
        for face in piece
            polygon = addpoint.(face)
            f = findfirst(f -> all(p -> onface(p, f), face), 1:6)
            if isnothing(f)
                # Membrane between the cell and the ECS
                push!(facets, polygon)
                push!(facetmarkers, interface_index(i, nregion, nregion))
                # ECS point: halfway to the neighboring cell, in the gap
                p = sum(face) / length(face) + gap / 2 * polygon_normal(face)
                all(lower .< p .< upper) && push!(ecsseeds, p)
            else
                push!(facets, polygon)
                push!(facetmarkers, outer_boundary_index(i, nregion))
                push!(boxpolygons[f], polygon)
            end
        end
    end

    # The faces of the box: the ECS parts are the faces with holes in the cells. The
    # outlines of the faces pass through the points of the cells on the edges of the box.
    for f = 1:6
        d = div(f + 1, 2)
        plane = isodd(f) ? lower[d] : upper[d]
        e1, e2 = filter(!=(d), 1:3)
        corners = [
            Base.setindex(Base.setindex(Base.setindex(lower, plane, d), x, e1), y, e2)
            for (x, y) in (
                (lower[e1], lower[e2]),
                (upper[e1], lower[e2]),
                (upper[e1], upper[e2]),
                (lower[e1], upper[e2]),
            )
        ]
        outline = addpoint.(corners)
        onedge = [
            k for (k, p) in enumerate(points) if onface(p, f) && (
                abs(p[e1] - lower[e1]) <= tol ||
                abs(p[e1] - upper[e1]) <= tol ||
                abs(p[e2] - lower[e2]) <= tol ||
                abs(p[e2] - upper[e2]) <= tol
            )
        ]
        center = (corners[1] + corners[3]) / 2
        angle(k) = atan(points[k][e2] - center[e2], points[k][e1] - center[e1])
        outline = sort(unique([outline; onedge]); by = angle)

        # Triangulate the face outside of the cells, without new points on the outline and
        # the cell polygons
        local2global = unique([outline; reduce(vcat, boxpolygons[f]; init = Int[])])
        global2local = Dict(k => i for (i, k) in enumerate(local2global))
        segments = Set{Tuple{Int,Int}}()
        for polygon in [[outline]; boxpolygons[f]], k in eachindex(polygon)
            a, b = global2local[polygon[k]],
            global2local[polygon[mod1(k + 1, length(polygon))]]
            push!(segments, minmax(a, b))
        end
        triin = TriangulateIO()
        triin.pointlist = reduce(hcat, [points[k][[e1, e2]] for k in local2global])
        triin.segmentlist = reduce(hcat, [collect(e) for e in segments])
        holes =
            [sum(points[polygon])[[e1, e2]] / length(polygon) for polygon in boxpolygons[f]]
        isempty(holes) || (triin.holelist = reduce(hcat, holes))
        triout = first(triangulate("pqYQ", triin))

        # New points inside the face
        for k = (length(local2global)+1):size(triout.pointlist, 2)
            q = triout.pointlist[:, k]
            p = Base.setindex(
                Base.setindex(Base.setindex(lower, plane, d), q[1], e1),
                q[2],
                e2,
            )
            push!(points, p)
            push!(local2global, length(points))
        end
        for t in eachcol(triout.trianglelist)
            push!(facets, local2global[t])
            push!(facetmarkers, outer_boundary_index(nregion, nregion))

            # ECS point: slightly inside the box from the triangle. The ECS in the box may
            # have several parts (e.g. pockets at the edges of the box), which all touch the
            # box faces.
            a, b, c = points[local2global[t]]
            area = norm(cross(b - a, c - a)) / 2
            inradius = 2 * area / (norm(b - a) + norm(c - b) + norm(a - c))
            inward = (isodd(f) ? 1 : -1) * SVector{3,T}(1:3 .== d)
            push!(ecsseeds, (a + b + c) / 3 + inradius / 1000 * inward)
        end
    end

    # Keep the ECS points outside the cells
    filter!(p -> !any(((_, piece),) -> inpolyhedron(p, piece), pieces), ecsseeds)
    isempty(ecsseeds) && error("No ECS found between the cells")
    append!(regions, ecsseeds)
    append!(regionmarkers, fill(nregion, length(ecsseeds)))

    (;
        points = reduce(hcat, points),
        facets,
        facetmarkers,
        regions = reduce(hcat, regions),
        regionmarkers,
    )
end

coefficients(setup::VoronoiSetup; kwargs...) =
    coefficients(SphereSetup(; ncell = setup.ncell, ecs = BoxECS()); kwargs...)
