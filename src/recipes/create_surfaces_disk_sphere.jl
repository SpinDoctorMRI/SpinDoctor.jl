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
            angles = 2π * collect(0:n-1) / n
            cell_points[j, i] = r * [cos.(angles'); sin.(angles')] .+ centers[:, j]
            cell_facets[j, i] = [(1:n-1)' n; (2:n)' 1]
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
