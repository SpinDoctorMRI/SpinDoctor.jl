"""
    create_surfaces(setup::CylinderSetup, cells)

Create surface triangulation of [inner and] outer cylinders [and ECS].

The ground surface is triangulated first, before the walls are "extruded" and
the top surface is copied from the ground surface.
"""
function create_surfaces(setup::CylinderSetup, cells)
    (; radii, centers) = cells
    (; ncell, r_range, height, include_in, in_ratio, ecs_shape, ecs_ratio) = setup
    (; rmin, rmax)  = r_range
    # Choose approximate cylinder side length
    nside = 30
    nside_min = 12

    include_ecs = ecs_shape != :no_ecs
    nboundary = (2 * include_in + 1 + include_ecs) * ncell + include_ecs
    rmean = (rmin + rmax) / 2

    # Create points on boundary circles
    create_nside(r) = max(nside_min, round(Int, nside * r / rmean))
    create_angles(n) = 2π * collect(0:n-1) / n
    create_circle(r, angles, center) = r * [cos.(angles'); sin.(angles')] .+ center
    create_circle_edges(n) = [(1:n-1)' n; (2:n)' 1]

    centers_columns = reshape([centers[:, i] for i = 1:ncell], (1, ncell))

    # Create points for IN compartments
    if include_in
        radii_in = in_ratio * radii
        nside_in = create_nside.(radii_in)
        angles_in = create_angles.(nside_in)
        circles_in = create_circle.(radii_in, angles_in, centers_columns)
        circle_edges_in = create_circle_edges.(nside_in)
        points_in = hcat(circles_in...)
        edges_in = hcat(circle_edges_in...)
        nprevious = 0
        for i = 1:ncell
            ncurrent = nprevious + nside_in[i]
            edges_in[:, nprevious+1:ncurrent] .+= nprevious
            nprevious = ncurrent
        end
    else
        # Create empty arrays
        nside_in = zeros(Int, 1, 0)
        points_in = zeros(2, 0)
        edges_in = zeros(Int, 2, 0)
    end


    # Create points of OUT compartments
    nside_out = create_nside.(radii)
    angles_out = create_angles.(nside_out)
    circles_out = create_circle.(radii, angles_out, centers_columns)
    circle_edges_out = create_circle_edges.(nside_out)
    points_out = hcat(circles_out...)
    edges_out = hcat(circle_edges_out...)
    nprevious = 0
    for i = 1:ncell
        ncurrent = nprevious + nside_out[i]
        edges_out[:, nprevious+1:ncurrent] .+= nprevious
        nprevious = ncurrent
    end

    if include_ecs
        radii_ecs = radii .+ ecs_ratio .* rmean
        nside_ecs = create_nside.(radii_ecs)
        angles_ecs = create_angles.(nside_ecs)
        circles_ecs = create_circle.(radii_ecs, angles_ecs, centers_columns)
        points_ecs = hcat(circles_ecs...)

        if ecs_shape == :box
            # Determine bounds of domain
            pmin = minimum(points_ecs, dims = 2)
            pmax = maximum(points_ecs, dims = 2)

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
            npoint_ecs = size(points_ecs, 2)
        elseif ecs_shape == :tight_wrap
            error("unimplemented")
            # ashape = alphashape(points_ecs)
            # edges_ecs = boundary(hull)
            # boundary_ecs = points_ecs[:, unique(edges_ecs)]
        end
        nedge_ecs = size(points_ecs, 2)
    else
        # Create empty arrays
        points_ecs = zeros(2, 0)
        edges_ecs = zeros(Int, 2, 0)
        nedge_ecs = zeros(Int, 1, 0)
    end

    # Create global collection of points and edges
    # Adjust numbering accordingly
    npoint_in = size(points_in, 2)
    npoint_out = size(points_out, 2)
    npoint_ecs = size(points_ecs, 2)
    points = [points_in points_out points_ecs]
    edges = [edges_in edges_out .+ npoint_in edges_ecs .+ npoint_in .+ npoint_out]
    npoint = size(points, 2)
    nedge = size(edges, 2)

    # Perform Delaynay triangulation of entire domain
    triin = TriangulateIO()
    triin.pointlist = points
    triin.segmentlist = edges
    triin.segmentmarkerlist = 1:nedge
    triout, = triangulate("pQ", triin)
    triangles = eachcol(triout.trianglelist)
    ntriangle = length(triangles)

    boundary_bounds = cumsum([0 nside_in nside_out nedge_ecs], dims = 2)

    ncell_in = ncell * include_in

    find_boundary(node) =
        findfirst(boundary_bounds[1:end-1] .+ 1 .≤ node .≤ boundary_bounds[2:end])

    # Create boundary numbers
    facetmarkers = zeros(Int, ntriangle)
    for (i, t) ∈ enumerate(triangles)
        b = find_boundary.(t)
        b_in = b .≤ ncell_in
        b_out = ncell_in .< b .≤ ncell_in + ncell
        if all(b_in)
            # Inside IN compartment
            facetmarkers[i] = ncell_in + include_ecs * ncell + b[1]
        elseif any(b_in)
            # Triangle touches IN compartment, but is not inside. It is thus the corresponding OUT compartment.
            facetmarkers[i] = ncell_in + include_ecs * ncell + b[.!b_in][1]
        elseif all(b_out) && b[1] == b[2] == b[3]
            # Triangle lies fully within the same OUT comparment
            facetmarkers[i] = ncell_in + include_ecs * ncell + b[1]
        else
            # Triangle lies in the ECS
            facetmarkers[i] = nboundary
        end
    end

    facets = hcat(triangles...)

    # Copy points, edges and markers two the two 3D planes top and bottom
    z = fill(height / 2, 1, npoint)
    points = [points points; -z z]

    function create_sidefacets(le, ri)
        [
            (le:ri-1)' ri (le+1:ri)'.+npoint le.+npoint
            (le:ri-1)'.+npoint ri.+npoint (le:ri-1)'.+npoint ri.+npoint
            (le+1:ri)' le (le+1:ri)' le
        ]
    end

    sidefacets =
        hcat(create_sidefacets.(boundary_bounds[1:end-1] .+ 1, boundary_bounds[2:end])...)

    sidefacetmarkers_in = 1:ncell_in
    sidefacetmarkers_in =
        vcat([repeat([m, m], nside_in[i]) for (i, m) ∈ enumerate(sidefacetmarkers_in)]...)

    sidefacetmarkers_out = (1+!include_ecs)ncell_in+1:(1+!include_ecs)ncell_in+ncell
    sidefacetmarkers_out =
        vcat([repeat([m, m], nside_out[i]) for (i, m) ∈ enumerate(sidefacetmarkers_out)]...)

    sidefacetmarkers_ecs = fill(nboundary, 2npoint_ecs)

    facets = [facets facets .+ npoint sidefacets]
    facetmarkers = [
        facetmarkers
        facetmarkers
        sidefacetmarkers_in
        sidefacetmarkers_out
        sidefacetmarkers_ecs
    ]

    if include_in
        regions_in = [centers; zeros(1, ncell)]
        regions_out = regions_in + [1; 0; 0] * (radii_in + radii) / 2
    else
        regions_in = zeros(3, 0)
        regions_out = [centers; zeros(1, ncell)]
    end
    if include_ecs
        xmin = argmin(points_ecs[1, :])
        region_ecs = [points_ecs[:, xmin]; 0] + [ecs_ratio * rmean / 2; 0; 0]
    else
        region_ecs = zeros(3, 0)
    end
    regions = [regions_in regions_out region_ecs]

    (; points, facets, facetmarkers, regions)
end
