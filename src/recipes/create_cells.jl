"""
    create_cells(setup)

Create geometrical configuration of cells. Return `nothing` or mathematical
description (radii, centers).
"""
function create_cells end

create_cells(::AbstractSetup) = nothing
create_cells(setup::ExtrusionSetup) = create_cells(setup.groundsetup)

function create_cells(setup::Union{DiskSetup,SphereSetup})
    (; ncell, rmin, rmax, dmin, dmax) = setup

    # Choose between spheres and cylinders
    if isa(setup, SphereSetup)
        ndim = 3
    elseif isa(setup, DiskSetup)
        ndim = 2
    end

    # Maximum number of points
    npoint_max = 100000 * ncell

    # Mean allowed radius
    rmean = (rmin + rmax) / 2

    # Cell centers (firsts is zero)
    centers = zeros(ndim, ncell)

    # Cell radii (first is rmean)
    radii = zeros(1, ncell)
    radii[1] = rmean

    # Initiate iterators
    icell = 1
    npoint = 1

    # Generate random cell centers until ncell cells are created
    while icell < ncell && npoint < npoint_max

        # Generate a random point using a uniform distribution with zero mean and
        # variance proportional to rmean
        if ndim == 3 # sphere
            point = (rand(ndim) .- 0.5) * rmean * max(10, ncell^(1 / 3))
        else # cylinder
            point = (rand(ndim) .- 0.5) * rmean * max(10.0, sqrt(ncell)) * 40.0
        end

        # Distance from point to the (already determined) cells
        # dist = sqrt.(sum(abs2, centers .- point, dims = 1)) .- radii
        # dist = dist[1:icell]
        dist = sqrt.(sum(abs2, centers[:, 1:icell] .- point, dims = 1)) .- radii[:, 1:icell]

        # Distance from point to nearest cell
        dist = minimum(dist)

        # Maximum allowed distance to nearest cell gives lowest allowed radius
        rmin1 = dist - dmax * rmean

        # Minimum allowed distance to nearest cell gives highest allowed radius
        rmax1 = dist - dmin * rmean

        # Radius has to be higher than both of the lowest allowed radii
        r_lower = max(rmin, rmin1)

        # Radius has to be lower than both of the highest allowed radii
        r_upper = min(rmax, rmax1)

        # Check if any radius is allowed
        if r_lower <= r_upper
            # Take midpoint between the bounds of all allowed radii
            r_use = (r_lower + r_upper) / 2

            # Update cells
            icell = icell + 1
            centers[:, icell] = point
            radii[icell] = r_use
        end

        # Update count for number of random centers generated
        npoint = npoint + 1
    end

    # Check that all cells were created
    icell == ncell || error("Did not find enough cell centers.")

    # Center cell collection
    pmin = minimum(centers .- radii, dims = 2)
    pmax = maximum(centers .+ radii, dims = 2)
    pmean = (pmin + pmax) / 2
    centers = centers .- pmean

    (; radii, centers)
end
