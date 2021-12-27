function get_volumes(setup)
    if isa(setup, CylinderSetup)
        dim = 2
        get_vol = r -> π * r^2 * setup.height
    else
        dim = 3
        get_vol = r -> 4 * π / 3 * r^3
    end

    # Get spherical volumes
    volumes = get_vol.(r * 1e6)

    # Subtract inner volumes (each volume is contained within the next)
    volumes[2:end] = volumes[2:end] - volumes[1:end-1]

    volumes
end
