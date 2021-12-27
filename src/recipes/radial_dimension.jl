"""
    radial_dimension(setup)

Get dimension of radial direction.
"""
function radial_dimension end

radial_dimension(::PlateSetup) = 1
radial_dimension(::CylinderSetup) = 2
radial_dimension(::SphereSetup) = 3
