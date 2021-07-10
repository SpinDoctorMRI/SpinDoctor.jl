"""
    elements, points = convexhull(points)

Compute the convex hull of a set of points `points` (of size `dim * npoint`).
Return a matrix of boundary elements `elements` (of size`dim * nelement`) and a restriction of the original points to the boundary (size `dim * npoint_keep`).
"""
function convexhull(points)

    # Sizes
    dim, npoint = size(points)

    # Compute Delaunay triangulation
    D = Int.(delaunay(points))

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
