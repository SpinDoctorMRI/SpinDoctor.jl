"""
    create_surfaces_sphere(cells, setup)

Create surface triangulation of [inner and] outer spheres [and ECS].
"""
function create_surfaces_sphere(cells, setup::Setup)

    # Extract parameters
    @unpack radii, centers = cells
    @unpack ncell, rmin, rmax, include_in, in_ratio, ecs_shape, ecs_ratio = setup

    include_ecs = ecs_shape != "no_ecs"

    rmean = (rmin + rmax) / 2

    create_npoint(radius) = round(Int, 200 * (radius / rmean)^2)
    if include_in
        npoint_in = create_npoint.(in_ratio .* radii)
        points_in = [
            centers[:, i] .+ in_ratio * radii[i] * create_fibonacci_sphere(npoint_in[i]) for i = 1:ncell
        ]
        facets_in = [convexhull(points_in[i])[1] for i = 1:ncell]
        regions_in = centers
        nfacet_in = size.(facets_in, 2)
    else
        npoint_in = zeros(Int, 0)
        nfacet_in = zeros(Int, 0)
        points_in = zeros(3, 0)
        facets_in = zeros(Int, 3, 0)
        regions_in = zeros(3, 0)
    end

    npoint_out = create_npoint.(radii)
    points_out =
        [centers[:, i] .+ radii[i] * create_fibonacci_sphere(npoint_out[i]) for i = 1:ncell]
    facets_out = [convexhull(points_out[i])[1] for i = 1:ncell]
    regions_out = @. centers + (1 + in_ratio) / 2 * radii * [1; 0; 0]
    nfacet_out = size.(facets_out, 2)

    if include_ecs
        if ecs_shape == "box"
            # Determine bounds of domain
            emin = minimum(hcat(points_out...), dims = 2)
            emax = maximum(hcat(points_out...), dims = 2)

            # Extend bounds by ECS gap
            @. emin = emin - ecs_ratio * rmean
            @. emax = emax + ecs_ratio * rmean

            # Define eight corners of box ECS and their corresponding edges
            points_ecs = [
                emin[1] emax[1] emax[1] emin[1] emin[1] emax[1] emax[1] emin[1]
                emin[2] emin[2] emax[2] emax[2] emin[2] emin[2] emax[2] emax[2]
                emin[3] emin[3] emin[3] emin[3] emax[3] emax[3] emax[3] emax[3]
            ]
            facets_ecs = [
                1 1 1 1 2 2 3 3 4 4 5 5
                2 3 2 6 3 7 4 8 1 5 6 7
                3 4 6 5 7 6 8 7 5 8 7 8
            ]
            regions_ecs = emin .+ ecs_ratio / 2 * rmean
            npoint_ecs = 8
            nfacet_ecs = 12
        elseif ecs_shape == "convex_hull"
            npoint_ecs = create_npoint.((1 + ecs_ratio) .* radii)
            points_ecs = [
                centers[:, i] .+
                (1 + ecs_ratio) * radii[i] * create_fibonacci_sphere(npoint_ecs[i]) for
                i = 1:ncell
            ]
            points_ecs = hcat(points_ecs...)

            facets_ecs, points_ecs = convexhull(points_ecs)
            npoint_ecs = size(points_ecs, 2)
            nfacet_ecs = size(facets_ecs, 2)
            regions_ecs = centers[:, 1] + (1 + ecs_ratio / 2) * radii[1] * [1; 0; 0]
        end
    else
        npoint_ecs = zeros(Int, 0)
        nfacet_ecs = zeros(Int, 0)
        points_ecs = zeros(3, 0)
        facets_ecs = zeros(Int, 3, 0)
        regions_ecs = zeros(3, 0)
    end

    npoint = cumsum(vcat(0, npoint_in..., npoint_out..., npoint_ecs))
    nfacet = cumsum(vcat(0, nfacet_in..., nfacet_out..., nfacet_ecs))
    points = hcat(points_in..., points_out..., points_ecs)
    facets = hcat(facets_in..., facets_out..., facets_ecs)
    facetmarkers = zeros(Int, size(facets, 2))
    for i = 1:length(npoint)-1
        inds = nfacet[i]+1:nfacet[i+1]
        facets[:, inds] .+= npoint[i]
        facetmarkers[inds] .= i
    end
    regions = [regions_in regions_out regions_ecs]

    (; points, facets, facetmarkers, regions)
end
