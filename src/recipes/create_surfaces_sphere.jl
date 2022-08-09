"""
    create_surfaces(setup::SphereSetup, cells)

Create surface triangulation of a set of multi-layered spherical cells
immersed in an ECS.

Returns a named tuple containing:

  - `points` (size `3 × npoint`),
  - `facets` (size `3 × nedge`),
  - `facetmarkers` (size `nfacet`),
  - `regions` (one point inside each region) (size `3 × nregion`)
"""
function create_surfaces(setup::SphereSetup, cells)
    (; radii, centers) = cells
    (; ncell, layersizes, rmin, rmax, nsidewall, ecs_shape, ecs_ratio) = setup

    # Minimum number of side facets
    nside_min = 20

    @assert issorted(layersizes)
    @assert all(s -> 0 < s ≤ 1, layersizes)

    nlayer = length(layersizes)
    include_ecs = ecs_shape != :no_ecs
    rmean = (rmin + rmax) / 2

    create_npoint(radius) = round(Int, nsidewall * (radius / rmean)^2)

    # Create empty arrays
    cell_points = [zeros(0, 0) for _ = 1:ncell, _ = 1:nlayer]
    cell_facets = [zeros(Int, 0, 0) for _ = 1:ncell, _ = 1:nlayer]
    cell_regions = [zeros(0) for _ = 1:ncell, _ = 1:nlayer]
    for i = 1:nlayer, j = 1:ncell
        r = layersizes[i] * radii[j]
        npoint = create_npoint(r)
        cell_points[j, i] = centers[:, j] .+ r .* create_fibonacci_sphere(npoint)
        cell_facets[j, i] = convexhull(cell_points[j, i])[1]
        if i == 1
            cell_regions[j, i] = centers[:, j]
        else
            cell_regions[j, i] =
                centers[:, j] + [1; 0; 0] * (layersizes[i-1] + layersizes[i]) / 2 * radii[j]
        end
    end

    if include_ecs
        r = radii .+ ecs_ratio .* rmean
        npoint = create_npoint.(r)
        ecs_points = reduce(
            hcat,
            centers[:, i] .+ r[i] * create_fibonacci_sphere(npoint[i]) for i = 1:ncell
        )
        xmin = argmin(ecs_points[1, :])
        ecs_regions = ecs_points[:, xmin] + [ecs_ratio * rmean / 2; 0; 0]
        if ecs_shape == :box
            # Determine bounds of domain
            emin = minimum(ecs_points; dims = 2)
            emax = maximum(ecs_points; dims = 2)

            # Define eight corners of box ECS and their corresponding edges
            ecs_points = [
                emin[1] emax[1] emax[1] emin[1] emin[1] emax[1] emax[1] emin[1]
                emin[2] emin[2] emax[2] emax[2] emin[2] emin[2] emax[2] emax[2]
                emin[3] emin[3] emin[3] emin[3] emax[3] emax[3] emax[3] emax[3]
            ]
            ecs_facets = [
                1 1 1 1 2 2 3 3 4 4 5 5
                2 3 2 6 3 7 4 8 1 5 6 7
                3 4 6 5 7 6 8 7 5 8 7 8
            ]
        elseif ecs_shape == :convex_hull
            ecs_facets, ecs_points = convexhull(ecs_points)
        elseif ecs_shape == :tight_wrap
            error("Not implemented")
        end
    else
        ecs_points = zeros(3, 0)
        ecs_facets = zeros(Int, 3, 0)
        ecs_regions = zeros(3, 0)
    end

    points = hcat(cell_points..., ecs_points)
    regions = hcat(cell_regions..., ecs_regions)
    nregion = size(regions, 2)

    # Take into account previous points in facets
    nprevious = 0
    for i = 1:nlayer, j = 1:ncell
        cell_facets[j, i] .+= nprevious
        nprevious += size(cell_points[j, i], 2)
    end
    ecs_facets .+= nprevious
    facets = hcat(cell_facets..., ecs_facets)

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
        elseif include_ecs
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
    append!(facetmarkers, fill((nregion - 1) * nregion ÷ 2 + nregion, size(ecs_facets, 2)))

    (; points, facets, facetmarkers, regions)
end
