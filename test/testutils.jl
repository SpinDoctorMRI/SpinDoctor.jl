"""
    get_setup(SetupType)

Get a preconfigured setup instance of type `SetupType`.
"""
function get_setup end

get_setup(S::Type{PlateSetup{T}}) where {T} = S(;
    name = "plates/someplates",
    width = 20.0,
    depth = 20.0,
    heights = fill(5.0, 4),
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
    ecs_shape = :convex_hull,
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
    ecs_shape = :convex_hull,
    ecs_ratio = 0.3,
    # refinement = 0.5,
)

get_setup(S::Type{NeuronSetup{T}}) where {T} = S(;
    name = "neuron",
    ecs_shape = :no_ecs,
    ecs_ratio = 0.3,
    # refinement = 0.5,
)

"""
    get_coeffs(setup)

Get preconfigured compartment coefficients for `setup`.
"""
function get_coeffs end

get_coeffs(setup::PlateSetup) = coefficients(
    setup;
    D = [0.002 * I(3) for _ = 1:length(setup.heights)],
    T₂ = [Inf for _ = 1:length(setup.heights)],
    ρ = [1.0 for _ = 1:length(setup.heights)],
    κ = (;
        interfaces = [1e-4 for _ = 1:length(setup.heights)-1],
        boundaries = [0.0 for _ = 1:length(setup.heights)],
    ),
    γ = 2.67513e-4,
)

get_coeffs(setup::CylinderSetup) = coefficients(
    setup;
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)

get_coeffs(setup::SphereSetup) = coefficients(
    setup;
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)

get_coeffs(setup::NeuronSetup) = coefficients(
    setup;
    D = (; neuron = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; neuron = Inf, ecs = Inf),
    ρ = (; n = 1.0, neuron = 1.0, ecs = 1.0),
    κ = (; neuron_ecs = 1e-4, neuron = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)
