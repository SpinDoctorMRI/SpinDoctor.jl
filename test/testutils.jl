"""
    get_setup(SetupType)

Get a preconfigured setup instance of type `SetupType`.
"""
function get_setup end

get_setup(S::Type{PlateSetup{T}}) where {T} = S(;
    name = "plates/someplates",
    width = 1.0,
    depth = 1.0,
    heights = 0:0.1:1.0,
    bend = 0.0,
    twist = 0.0,
    # refinement = 0.5,
)

get_setup(S::Type{CylinderSetup{T}}) where {T} = S(;
    name = "cylinders/somecylinders",
    ncell = 1,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    height = 1.0,
    bend = 0.0,
    twist = 0.0,
    include_in = true,
    in_ratio = 0.6,
    ecs_shape = "convex_hull",
    ecs_ratio = 0.5,
    # refinement = 0.5,
)

get_setup(S::Type{SphereSetup{T}}) where {T} = S(;
    name = "spheres/somespheres",
    ncell = 4,
    rmin = 3.0,
    rmax = 5.0,
    dmin = 0.2,
    dmax = 0.3,
    include_in = false,
    in_ratio = 0.7,
    ecs_shape = "convex_hull",
    ecs_ratio = 0.3,
    # refinement = 0.5,
)

get_setup(S::Type{NeuronSetup{T}}) where {T} = S(;
    name,
    ecs_shape = "no_ecs",
    ecs_ratio = 0.3,
    # refinement = 0.5,
)

"""
    setup_coeffs(setup)

Get preconfigured compartment coefficients for `setup`.
"""
setup_coeffs(setup) = coefficients(
    setup;
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)

