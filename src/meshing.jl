"""
    create_mesh(surfaces, refinement)

Create mesh from surface geometry. A refinement is applied if `refinement < Inf`.
For 3D, TetGen is called. For 2D, Triangle is called.
"""
function create_mesh(surfaces, refinement)
    (; points, facets, facetmarkers, regions) = surfaces

    dim = size(points, 1)
    nregion = size(regions, 2)

    if dim == 2
        regionlist = [regions; (1:nregion)'; fill(refinement, 1, nregion)]

        triin = TriangulateIO()
        triin.pointlist = points
        triin.segmentlist = facets
        triin.segmentmarkerlist = facetmarkers
        triin.regionlist = regionlist
        triout = first(triangulate("pqcaAQ", triin))

        points = triout.pointlist
        facets = Int.(triout.segmentlist)
        facetmarkers = Int.(triout.segmentmarkerlist)
        elements = Int.(triout.trianglelist)
        elementmarkers = round.(Int, triout.triangleattributelist[:])
    elseif dim == 3
        input = RawTetGenIO{Cdouble}()
        input.pointlist = points

        facetlist!(input, facets)

        input.facetmarkerlist = facetmarkers

        input.regionlist = [
            regions
            collect(1:nregion)'
            fill(0.1, 1, nregion)
        ]

        if isinf(refinement)
            tetgen = tetrahedralize(input, "pqAQ")
        else
            tetgen = tetrahedralize(input, "pqAQa$(refinement)")
        end

        points = tetgen.pointlist
        facets = Int.(tetgen.trifacelist)
        facetmarkers = Int.(tetgen.trifacemarkerlist)
        elements = Int.(tetgen.tetrahedronlist)
        elementmarkers = Int.(tetgen.tetrahedronattributelist[:])
    end

    (; points, facets, facetmarkers, elements, elementmarkers)
end

"""
    convexhull2(points)

Convex hull of 2D points (gift wrapping).
https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
"""
function convexhull2(p)
    p = copy(p)

    function swap!(points, i, j)
        tmp = points[:, i]
        points[:, i] = points[:, j]
        points[:, j] = tmp
    end

    isclock(u, v) = u[2]v[1] - u[1]v[2] ≥ 0

    npoint = size(p, 2)

    # First point with lowest x-coordinate
    swap!(p, 1, argmin(p[1, :]))

    # First u (turn clockwise from here)
    u = [0, 1]

    # Iterate through points
    i = 1
    for i = 1:(npoint-1)
        done = true
        jnext = i
        for j = (i+1):npoint
            if isclock(u, p[:, j] - p[:, i])
                u = p[:, j] - p[:, i]
                jnext = j
                done = false
            end
        end
        swap!(p, i + 1, jnext)
        u = p[:, i] - p[:, i+1]

        if done || (i ≥ 2 && isclock(u, p[:, i] - p[:, 1]))
            p = p[:, 1:i]
            break
        end
    end

    n = size(p, 2)
    e = [(1:n)'; (2:n)' 1]

    e, p
end

"""
    convexhull3(points)

Convex hull of 3D points. Based on
https://github.com/rileyborgard/convex-hull-3d/blob/master/convexhull3d.cpp
"""
function convexhull3(points)
    n = size(points, 2)
    @assert n ≥ 4

    points = SVector{3}.(eachcol(points))
    ϵ = 1e-10

    tet = [1]
    for i = 2:n
        length(tet) < 4 || break
        a = length(tet) == 1 && norm(points[tet[1]] - points[i]) > ϵ
        b =
            length(tet) == 2 &&
            norm((points[tet[2]] - points[tet[1]]) × (points[i] - points[tet[1]])) > ϵ
        c =
            length(tet) == 3 &&
            abs(
                (points[i] - points[tet[1]]) ⋅
                cross(points[tet[2]] - points[tet[1]], points[tet[3]] - points[tet[1]]),
            ) > ϵ
        if a || b || c
            push!(tet, i)
        end
    end
    @assert length(tet) == 4

    pointsfirst = points[tet]
    prepend!(deleteat!(points, tet), pointsfirst)

    faces = SVector{3,Int}[]

    dead = trues(n, n)
    function addface(a, b, c)
        push!(faces, SVector{3,Int}(a, b, c))
        dead[a, b] = dead[b, c] = dead[c, a] = false
    end

    addface(1, 2, 3)
    addface(1, 3, 2)

    for i = 4:n
        faces2 = SVector{3,Int}[]
        for f ∈ faces
            a, b, c = f[1], f[2], f[3]
            n = (points[b] - points[a]) × (points[c] - points[a])
            n /= norm(n)
            if (points[i] - points[a]) ⋅ n > ϵ
                dead[a, b] = dead[b, c] = dead[c, a] = true
            else
                push!(faces2, f)
            end
        end
        empty!(faces)
        for f ∈ faces2
            a, b, c = f[1], f[2], f[3]
            dead[b, a] && addface(b, a, i)
            dead[c, b] && addface(c, b, i)
            dead[a, c] && addface(a, c, i)
        end
        append!(faces, faces2)
    end

    faces = reduce(hcat, faces)
    points = reduce(hcat, points)

    faces, points
end

"""
    elements, points = convexhull(points)

Compute the convex hull of a set of points `points` (of size `dim * npoint`). Return a
matrix of boundary elements `elements` (of size`dim * nelement`) and a restriction of the
original points to the boundary (size `dim * npoint_keep`).
"""
function convexhull(points)
    dim = size(points, 1)
    if dim == 2
        e, p = convexhull2(points)
    elseif dim == 3
        e, p = convexhull3(points)
    else
        error("convexhull only implemented for 2D and 3D")
    end
    e, p
end

# # 2D
# a = rand(2, 100)
# e, p = convexhull(a)
# pl = scatter(a[1, :], a[2, :])
# scatter!(pl, p[1, :], p[2, :], color = :red)
# for i = 1:size(e, 2)
#     plot!(pl, p[1, e[:, i]], p[2, e[:, i]])
#     display(pl)
#     sleep(0.5)
# end

# # 3D
# points = randn(3, 1000)
# e, p = convexhull3(points)
# poly(p', e'; color = p[3, :], transparency = true, strokewidth = 1, shading = false)
# scatter!(points)

# Deform domain by bending and twisting.
function deform_domain!(points, bend, twist)
    height = maximum(points[3, :]) - minimum(points[3, :])
    width = maximum(points[1, :]) - minimum(points[1, :])

    # Ben points in x-direction as a percentage of height.
    zcenter = (maximum(points[3, :]) + minimum(points[3, :])) / 2
    @. points[1, :] += bend * 30.0 * width * ((points[3, :] - zcenter) / height)^2.0

    # Twist angles by height
    thvec = points[[3], :] / height * twist

    # Twist around center
    center = (maximum(points[1:2, :], dims = 2) .+ minimum(points[1:2, :], dims = 2)) ./ 2
    points[1:2, :] .-= center

    points[1:2, :] .= [
        cos.(thvec) .* points[[1], :] .- sin.(thvec) .* points[[2], :]
        sin.(thvec) .* points[[1], :] .+ cos.(thvec) .* points[[2], :]
    ]

    # points[1, :] = points[1, :] - rmax
    points[1:2, :] .+= center

    nothing
end
