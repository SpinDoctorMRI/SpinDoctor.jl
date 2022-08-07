"""
    create_surfaces(setup::DiskSetup, cells)

Create surface triangulation of [inner and] outer cylinders [and ECS].

The surface is triangulated.
"""
function create_surfaces(setup::DiskSetup, cells)
    (; radii, centers) = cells
    (; ncell, nsidewall, rmin, rmax, layersizes, ecs_shape, ecs_ratio, refinement) = setup

    # Choose approximate cylinder side length
    nside_min = 8

    @assert issorted(layersizes)
    @assert all(s -> 0 < s ≤ 1, layersizes)

    nlayer = length(layersizes)
    include_ecs = ecs_shape != :no_ecs
    rmean = (rmin + rmax) / 2

    # Create points on boundary circles
    create_nside(r) = max(nside_min, round(Int, nsidewall * r / rmean))
    create_angles(n) = 2π * collect(0:n-1) / n
    create_circle(r, angles, center) = r * [cos.(angles'); sin.(angles')] .+ center
    create_circle_edges(n) = [(1:n-1)' n; (2:n)' 1]

    centers_columns = reshape([centers[:, i] for i = 1:ncell], (1, ncell))

    # Create points for IN compartments
    # Create empty arrays
    points_in = [zeros(0, 0) for _ = 1:ncell, _ = 1:nlayer]
    edges_in = [zeros(Int, 0, 0) for _ = 1:ncell, _ = 1:nlayer]
    regions_in = [zeros(0) for _ = 1:ncell, _ = 1:nlayer]
    for i = 1:nlayer, j = 1:ncell
        r = layersizes[i] * radii[j]
        nside = create_nside(r)
        angles = create_angles(nside)
        points_in[j, i] = create_circle(r, angles, centers[:, j])
        edges_in[j, i] = create_circle_edges(nside)
        if i == 1
            regions_in[j, i] = centers[:, j]
        else
            regions_in[j, i] =
                centers[:, j] + [1; 0] * (layersizes[i-1] + layersizes[i]) / 2 * radii[j]
        end
    end

    if include_ecs
        radii_ecs = radii .+ ecs_ratio .* rmean
        nside_ecs = create_nside.(radii_ecs)
        angles_ecs = create_angles.(nside_ecs)
        circles_ecs = create_circle.(radii_ecs, angles_ecs, centers_columns)
        points_ecs = hcat(circles_ecs...)

        if ecs_shape == :box
            # Determine bounds of domain
            pmin = minimum(points_ecs; dims = 2)
            pmax = maximum(points_ecs; dims = 2)

            # Define four corners of box ECS and their corresponding edges
            points_ecs = [
                pmin[1] pmax[1] pmax[1] pmin[1]
                pmin[2] pmin[2] pmax[2] pmax[2]
            ]
            edges_ecs = [
                1 2 3 4
                2 3 4 1
            ]
        elseif ecs_shape == :convex_hull
            # Extract points defining convex hull of ECS
            edges_ecs, points_ecs = convexhull(points_ecs)
        elseif ecs_shape == :tight_wrap
            error("unimplemented")
            # ashape = alphashape(points_ecs)
            # edges_ecs = boundary(hull)
            # boundary_ecs = points_ecs[:, unique(edges_ecs)]
        end
        xmin = argmin(points_ecs[1, :])
        region_ecs = points_ecs[:, xmin] + [ecs_ratio * rmean / 2; 0]
    else
        # Create empty arrays
        points_ecs = zeros(2, 0)
        edges_ecs = zeros(Int, 2, 0)
        region_ecs = zeros(2, 0)
    end

    points = hcat(points_in..., points_ecs)
    regions = hcat(regions_in..., region_ecs)
    nregion = size(regions, 2)

    # Take into account previous points in edges
    nprevious = 0
    for i = 1:nlayer, j = 1:ncell
        edges_in[j, i] .+= nprevious
        nprevious += size(points_in[j, i], 2)
    end
    edges_ecs .+= nprevious
    edges = hcat(edges_in..., edges_ecs)

    # The boundaries are ordered as follows:
    #   Number    Boundary
    #      1       Γ_1_2 (we skip Γ_2_1)
    #      2       Γ_1_3 (we skip Γ_3_1)
    #      ⋮         ⋮
    #  (n-1)n/2    Γ_n_n
    # (n-1)n/2+1    Γ_1
    #                ⋮
    # (n-1)n/2+n    Γ_n
    edgenumber(i, j) = (i - 1) * (2 * nregion - i) ÷ 2 + (j - i)
    edgemarkers = zeros(Int, 0)
    k = 1
    for i = 1:nlayer, j = 1:ncell
        if i < nlayer
            # (i, j) is linked to (i + 1, j)
            m = edgenumber(k, k + ncell)
        elseif include_ecs
            # The outer layer is linked to the ECS (last region)
            m = edgenumber(k, nregion)
        else
            # The outer region k has an outer boundary Γ_k
            m = (nregion - 1) * nregion ÷ 2 + k
        end
        append!(edgemarkers, fill(m, size(edges_in[j, i], 2)))
        k += 1
    end

    # ECS is last region, with outer boundary Γ_nregion
    append!(edgemarkers, fill((nregion - 1) * nregion ÷ 2 + nregion, size(edges_ecs, 2)))

    (; points, facets = edges, facetmarkers = edgemarkers, regions)
end
