function create_surfaces(setup::PlateSetup{T}, _) where {T}
    (; widths, depth) = setup

    n = length(widths)
    w = sum(widths) / 2
    x = cumsum([-w; widths])
    y = depth / 2

    points = zeros(T, 2, 2(n + 1))
    for i = 1:n+1
        points[:, 1+2(i-1):2i] = [
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
