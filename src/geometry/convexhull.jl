"""
    elements, points = convexhull(points)

Compute the convex hull of a set of points `points` (of size `dim * npoint`). Return a
matrix of boundary elements `elements` (of size`dim * nelement`) and a restriction of the
original points to the boundary (size `dim * npoint_keep`).
"""
function convexhull(points)

    # Sizes
    dim, npoint = size(points)

    # Compute Delaunay triangulation
    if dim == 2
        return convexhull2(points)
    else
        # D = Int.(delaunay3(points))
        D = Int.(delaunay(points))
    end

    # Sort each column of d
    for i = 1:size(D, 2)
        sort!(@view D[:, i])
    end

    # Identify boundary elements. In 3D, those are facets not shared by two elements, in 2D
    # those are edges not shared by two facets
    if dim == 2
        combs = [[1, 2], [1, 3], [2, 3]]
    else
        combs = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
    end
    A = [SVector{dim}(D[comb, i]) for i = 1:size(D, 2) for comb ∈ combs]
    sort!(A)

    # Indices of unique elements
    ikeep = fill(false, length(A))

    # Perform lookup for the first element
    lookup = true

    # Loop through facets or edges. Each one occurs exactly once or twice, and they are
    # ordered
    for i ∈ 1:length(A)-1
        if lookup
            # A[i] is different than A[i-1]
            if A[i] == A[i+1]
                # A[i] and A[i+1] is a double element, and should not be kept. Skip next
                # iteration
                lookup = false
            else
                # A[i] is unique, and is thus a boundary element
                ikeep[i] = true
            end
        else
            # A[i] is a copy of A[i-1], but A[i+1] is a new element. Perform next iteration
            lookup = true
        end
    end

    # Keep the last element if it is different from the one before
    ikeep[end] = lookup

    # Keep unique elements
    A = A[ikeep]

    # Order edges to form a connected wrapping in 2D case
    if dim == 2
        for i = 1:length(A)-1
            # Find next edge from endpoint of previous edge
            j = i + findfirst(e -> A[i][2] ∈ e, @view A[i+1:end])

            # Reorder
            tmp = A[i+1]
            A[i+1] = A[j]
            A[j] = tmp

            # Orient next edge correctly
            if A[i+1][2] == A[i][2]
                A[i+1] = SVector{2}(A[i+1][[2, 1]])
            end
        end
    end

    E = mapreduce(Vector, hcat, A)

    # Important: `unique` preserves order of edges in 2D case
    ikeep = unique(E)

    # Inverse point map
    ikeep_inv = fill(0, npoint)
    for i = 1:length(ikeep)
        ikeep_inv[ikeep[i]] = i
    end

    # Renumber points from 1 to nkeep
    boundary_elements = ikeep_inv[E]
    boundary_points = points[:, ikeep]

    boundary_elements, boundary_points
end

# a = rand(2, 100)
# e, p = convexhull(a)

# pl = scatter(a[1, :], a[2, :])

# scatter!(pl, p[1, :], p[2, :], color = :red)

# for i = 1:size(e, 2)
#     plot!(pl, p[1, e[:, i]], p[2, e[:, i]])
#     display(pl)
#     sleep(0.5)
# end


"""
    delaunay3(points)

3D Delaunay triangulation from Tetgen.
"""
function delaunay3(points)
    input = RawTetGenIO{Cdouble}()
    input.pointlist = points
    tetgen = tetrahedralize(input, "cQ")

    tetgen.tetrahedronlist
end


"""
    convexhull2(points)

Convex hull of 2D points (gift wrapping).
https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
"""
function convexhull2(p)
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
    for i = 1:npoint-1
        done = true
        jnext = i
        for j = i+1:npoint
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
