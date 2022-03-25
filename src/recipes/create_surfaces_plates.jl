"""
    create_surfaces(setup::PlateSetup, cells)

Create surface triangulation of a collection of stacked plates.
"""
function create_surfaces(setup::PlateSetup{T}, _) where {T}
    (; width, depth, heights) = setup

    n = length(heights)
    h = sum(heights) / 2
    z = cumsum([-h; heights])
    x = width / 2
    y = depth / 2

    points = zeros(T, 3, 4(n + 1))
    for i = 1:(n+1)
        points[:, (1+4(i-1)):(4i)] = [
            -x x x -x
            -y -y y y
            z[i] z[i] z[i] z[i]
        ]
    end

    facets = fill(0, 3, 2 + 10n)
    facets[:, 1:2] = [
        1 1
        2 3
        3 4
    ]
    for i = 1:n
        inds = 2 + 10(i - 1) .+ (1:10)
        facets[:, inds] .=
            4(i - 1) .+ [
                1 1 2 2 3 3 4 4 6 6
                2 6 3 7 4 8 1 5 7 8
                6 5 7 6 8 7 5 8 8 5
            ]
    end

    facetmarkers = fill(0, 2 + 10n)

    # Bottom and first side walls
    facetmarkers[1:10] .= n

    # Interfaces and side walls after interface
    for i = 1:(n-1)
        facetmarkers[10i+1] = i
        facetmarkers[10i+2] = i
        facetmarkers[(10i+3):(10(i+1))] .= n + i
    end

    # Top
    facetmarkers[(end-1):end] .= 2n - 1

    regions = zeros(T, 3, n)
    regions[3, :] .= (z[1:(end-1)] .+ z[2:end]) ./ 2


    (; points, facets, facetmarkers, regions)
end
