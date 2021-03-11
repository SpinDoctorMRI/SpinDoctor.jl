""" Create surface triangulation of spheres or cylinders, with or without nuclei and ECS. """
function create_surface_triangulation(cellsetup::CellSetup, cells, domain)
    @unpack shape = cellsetup
    if shape == "cylinder"
        tri = create_surface_triangulation_cylinder(cellsetup, cells, domain)
    elseif shape == "sphere"
        tri = create_surface_triangulation_sphere(cellsetup, cells, domain)
    end
    tri
end


"""
    create_surface_triangulation(cellsetup, cells, domain)

Create surface triangulation of [inner and] outer cylinders [and ECS].

The ground surface is triangulated first, before the walls are "extruded" and
the top surface is copied from the ground surface.
"""
function create_surface_triangulation_cylinder(cellsetup::CellSetup, cells, domain)
    @unpack ncell, rmin, rmax, height, include_in, in_ratio, include_ecs, ecs_shape, ecs_ratio = cellsetup
    @unpack radii, centers = cells
    @unpack compartments, boundaries, ncompartment, nboundary = domain

    # Choose approximate cylinder side length
    nside = 30
    nside_min = 12
    rmean = (rmin + rmax) / 2

    # Create points on boundary circles
    create_nside(r) = max(nside_min, round(Int, nside * r / rmean))
    create_angles(n) = 2π * (0:n-1) / n
    create_circle(r, angles, center) = r * [cos.(angles'); sin.(angles')] .+ center
    create_circle_edges(n) = [collect(1:n-1)' n; collect(2:n)' 1]

    centers_columns = reshape([centers[:, i] for i=1:ncell], (1, ncell))

    # Create points of OUT compartments
    nside_out = create_nside.(radii)
    angles_out = create_angles.(nside_out)
    circles_out = create_circle.(radii, angles_out, centers_columns)
    circle_edges_out = create_circle_edges.(nside_out)

    points_out = hcat(circles_out...)'
    edges_out = hcat(circle_edges_out...)'

    nprevious = 0
    for i = 1:ncell
        ncurrent = nprevious + nside_out[i]
        edges_out[nprevious+1:ncurrent, :] .+= nprevious
        nprevious = ncurrent
    end

    # Create points for IN compartments
    if include_in
        radii_in = in_ratio * radii
        nside_in = create_nside.(radii_in)
        angles_in = create_angles.(nside_in)
        circles_in = create_circle.(radii_in, angles_in, centers_columns)
        circle_edges_in = create_circle_edges.(nside_in)

        points_in = hcat(circles_in...)'
        edges_in = hcat(circle_edges_in...)'
        nprevious = 0
        for i = 1:ncell
            ncurrent = nprevious + nside_in[i]
            edges_in[nprevious+1:ncurrent, :] .+= nprevious
            nprevious = ncurrent
        end
    else
        # Create empty arrays
        nside_in = zeros(Int, 1, 0)
        points_in = zeros(0, 2)
        edges_in = zeros(Int, 0, 2)
    end

    if include_ecs
        radii_ecs = radii .+ ecs_ratio * rmean
        nside_ecs = create_nside.(radii_ecs)
        angles_ecs = create_angles.(nside_ecs)
        circles_ecs = create_circle.(radii_ecs, angles_ecs, centers_columns)
        points_ecs = hcat(circles_ecs...)'

        if ecs_shape == "box"
            # Determine bounds of domain
            emin = minimum(hcat(circles_out...), dims=2)
            emax = maximum(hcat(circles_out...), dims=2)

            # Extend bounds by ECS gap
            @. emin = emin - ecs_ratio * rmean
            @. emax = emax + ecs_ratio * rmean

            # Define four corners of box ECS and their corresponding edges
            points_ecs = [
                emin[1] emin[2]
                emax[1] emin[2]
                emax[1] emax[2]
                emin[1] emax[2]
            ]
            edges_ecs = [1 2; 2 3; 3 4; 4 1]
            nedge_ecs = 4

        elseif ecs_shape == "convexhull"
            # Extract points defining convex hull of ECS
            lib = DefaultLibrary{Float64}(Optimizer)
            polyhed = polyhedron(vrep(points_ecs), lib)
            removevredundancy!(polyhed)
            points_ecs = polyhed.vrep.V
            nedge_ecs = size(points_ecs, 1)
            edges_ecs = [1:nedge_ecs-1 2:nedge_ecs; nedge_ecs 1]
        elseif ecs_shape == "tight_wrap"
            error("unimplemented")
            # ashape = alphashape(points_ecs)
            # edges_ecs = boundary(hull)
            # boundary_ecs = points_ecs[:, unique(edges_ecs)]
        end
        nside_ecs = size(points_ecs, 1)
    else
        # Create empty arrays
        nside_ecs = zeros(Int, 1, 0)
        points_ecs = zeros(0, 2)
        edges_ecs = zeros(Int, 0, 2)
    end

    # Create global collection of points and edges
    # Adjust numbering accordingly
    npoint_in = size(points_in, 1)
    npoint_out = size(points_out, 1)
    npoint_ecs = size(points_ecs, 1)

    nedge_in = size(edges_in, 1)
    nedge_out = size(edges_out, 1)
    nedge_ecs = size(edges_ecs, 1)

    points = [
        points_in
        points_out
        points_ecs
    ]
    # edges = [
    #     edges_in
    #     edges_out .+ nedge_in
    #     edges_ecs .+ nedge_in .+ nedge_out
    # ]
    edges = [
        edges_in
        edges_out .+ npoint_in
        edges_ecs .+ npoint_in .+ npoint_out
    ]
    npoint = size(points, 1)
    nedge = size(edges, 1)

    # Perform Delaynay triangulation of entire domain
    triangles = constrained_triangulation(points, collect(1:npoint), edges, repeat([true], nedge))
    ntriangle = length(triangles)

    boundary_bounds = cumsum([0 nside_in nside_out nside_ecs], dims=2)

    ncell_in = ncell * include_in

    find_boundary(node) = findfirst(boundary_bounds[1:end-1] .< node .≤ boundary_bounds[2:end])

    # Create boundary numbers
    trianglemarkers = zeros(Int, ntriangle)
    for (i, t) ∈ enumerate(triangles)
        b = find_boundary.(t)
        b_in = b .≤ ncell_in
        b_out = ncell_in .< b .≤ ncell_in + ncell
        if all(b_in)
            # Inside IN compartment
            trianglemarkers[i] = b[1]
        elseif any(b_in)
            # Triangle touches IN compartment, but is not inside. It is thus the corresponding OUT compartment.
            trianglemarkers[i] = b[.!b_in][1]
        elseif all(b_out) && b[1] == b[2] == b[3]
            # Triangle lies fully within the same OUT comparment
            trianglemarkers[i] = b[1]
        else
            # Triangle lies in the ECS
            trianglemarkers[i] = nboundary
        end
    end

    facets = hcat(triangles...)'

    # Copy points, edges and markers two the two 3D planes top and bottom
    z = repeat([height / 2], npoint)
    points3d = [points -z; points z]

    function create_sidefacets(le, ri)
        [
                le:ri-1       (le:ri-1).+npoint  (le+1:ri)
                  ri              ri+npoint         le
            (le+1:ri).+npoint (le:ri-1).+npoint  (le+1:ri)
               le+npoint          ri+npoint         le
        ]
    end

    sidefacets = vcat(create_sidefacets.(boundary_bounds[1:end-1] .+ 1, boundary_bounds[2:end])...)

    sidefacetmarkers_in = ncell_in+ncell+1:2ncell_in+ncell
    sidefacetmarkers_in = vcat([repeat([m, m], nside_in[i]) for (i, m) ∈ enumerate(sidefacetmarkers_in)]...)

    sidefacetmarkers_out = 2ncell_in+ncell+1:2ncell_in+2ncell
    sidefacetmarkers_out = vcat([repeat([m, m], nside_out[i]) for (i, m) ∈ enumerate(sidefacetmarkers_out)]...)

    sidefacetmarkers_ecs = repeat([nboundary], 2npoint_ecs)

    facets3d = [facets; facets .+ npoint; sidefacets]
    facetmarkers = [trianglemarkers; trianglemarkers; sidefacetmarkers_in; sidefacetmarkers_out; sidefacetmarkers_ecs]

    if include_in
        regions_in = [1 0; 0 1; 0 0] * centers
        regions_out = regions_in + [1; 0; 0] * (radii_in + radii) / 2
    else
        regions_in = zeros(3, 0)
        regions_out = [1 0; 0 1; 0 0] * centers
    end
    if include_ecs
        xmin = argmin(points_ecs[:, 1])
        region_ecs = [points_ecs[xmin, :]; 0] + [ecs_ratio*rmean/2; 0; 0]
    else
        region_ecs = zeros(3, 0)
    end
    regions = [regions_in'; regions_out'; region_ecs']

    points = points3d
    facets = facets3d

    points, facets, regions = points', facets', regions'

    (; points, facets, facetmarkers, regions)
end


"""
    create_surface_triangulation_sphere(cellsetup, cells, domain)

Create surface triangulation of [inner and] outer spheres [and ECS].
"""
function create_surface_triangulation_sphere(cellsetup::CellSetup, cells, domain)

    # Extract parameters
    @unpack ncell, rmin, rmax, height, include_in, in_ratio, include_ecs, ecs_shape, ecs_ratio = cellsetup
    @unpack radii, centers = cells
    @unpack compartments, boundaries, ncompartment, nboundary = domain

    rmean = (rmin + rmax) / 2

    create_npoint(radius) = round(Int, 200 * (radius / rmean)^2)
    if include_in
        npoint_in = create_npoint.(in_ratio .* radii)
        points_in = [centers[:, i] .+ in_ratio * radii[i] * create_fibonacci_sphere(npoint_in[i]) for i = 1:ncell]
        facets_in = [hcat(chull(collect(points_in[i]')).simplices...) for i = 1:ncell]
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
    points_out = [centers[:, i] .+ radii[i] * create_fibonacci_sphere(npoint_out[i]) for i = 1:ncell]
    facets_out = [hcat(chull(collect(points_out[i]')).simplices...) for i = 1:ncell]
    regions_out = @. centers + (1 + in_ratio) / 2 * radii * [1; 0; 0]
    nfacet_out = size.(facets_out, 2)
    
    if include_ecs
        if ecs_shape == "box"
            # Determine bounds of domain
            emin = minimum(hcat(points_out...), dims=2)
            emax = maximum(hcat(points_out...), dims=2)

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
            regions_ecs =  emin .+ ecs_ratio / 2 * rmean
            npoint_ecs = 8
            nfacet_ecs = 12
        elseif ecs_shape == "convexhull"
            npoint_ecs = create_npoint.((1 + ecs_ratio) .* radii)
            points_ecs = [centers[:, i] .+ (1 + ecs_ratio) * radii[i] * create_fibonacci_sphere(npoint_ecs[i]) for i = 1:ncell]
            points_ecs = hcat(points_ecs...)

            ecs_hull = chull(collect(points_ecs'))
            points_ecs = points_ecs[:, ecs_hull.vertices]
            npoint_ecs = size(points_ecs, 2)
            facets_ecs = hcat(ecs_hull.simplices...)
            [facets_ecs[i, j] = findfirst(facets_ecs[i, j] .== ecs_hull.vertices) for j = 1:length(ecs_hull.simplices) for i = 1:3]
            nfacet_ecs = size(facets_ecs, 2)
            regions_ecs =  centers[:, 1] + (1 + ecs_ratio / 2) * radii[1] * [1; 0; 0]
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
