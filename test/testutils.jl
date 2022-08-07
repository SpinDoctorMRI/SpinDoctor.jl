"""
    get_setup(SetupType)

Get a preconfigured setup instance of type `SetupType`.
"""
function get_setup end

get_setup(S::Type{PlateSetup{T}}) where {T} = S(;
    depth = 20.0,
    widths = fill(5.0, 4),
    # refinement = 0.5,
)

get_setup(S::Type{DiskSetup{T}}) where {T} = S(;
    ncell = 1,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    layersizes = [0.2, 0.7, 1.0],
    ecs_shape = :convex_hull,
    ecs_ratio = 0.5,
    # refinement = 0.5,
)

get_setup(S::Type{SlabSetup{T}}) where {T} = S(;
    depth = 20.0,
    widths = fill(5.0, 4),
    height = 20.0,
    bend = 0.0,
    twist = 0.0,
    # refinement = 0.5,
)

get_setup(S::Type{CylinderSetup{T}}) where {T} = S(;
    ncell = 1,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    layersizes = [0.2, 0.7, 1.0],
    ecs_shape = :convex_hull,
    ecs_ratio = 0.5,
    height = 5,
    bend = 0,
    twist = 0,
    # refinement = 0.5,
)

get_setup(S::Type{SphereSetup{T}}) where {T} = S(;
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
    ecs_shape = :no_ecs,
    # refinement = 0.5,
)

"""
    get_coeffs(setup)

Get preconfigured compartment coefficients for `setup`.
"""
function get_coeffs end

function get_coeffs(setup::PlateSetup)
    nlayer = length(setup.widths)
    coefficients(
        setup;
        D = [0.002 * I(2) for _ = 1:nlayer],
        T₂ = fill(Inf, nlayer),
        ρ = fill(1.0, nlayer),
        κ = (; interfaces = fill(1e-4, nlayer), boundaries = fill(0.0, nlayer)),
        γ = 2.67513e-4,
    )
end

function get_coeffs(setup::SlabSetup)
    nlayer = length(setup.groundsetup.widths)
    coefficients(
        setup;
        D = [0.002 * I(3) for _ = 1:nlayer],
        T₂ = fill(Inf, nlayer),
        ρ = fill(1.0, nlayer),
        κ = (; interfaces = fill(1e-4, nlayer), boundaries = fill(0.0, nlayer)),
        γ = 2.67513e-4,
    )
end

function get_coeffs(setup::DiskSetup)
    nlayer = length(setup.layersizes)
    coefficients(
        setup;
        D = (; cell = [0.002 * I(2) for _ = 1:nlayer], ecs = 0.002 * I(2)),
        T₂ = (; cell = fill(Inf, nlayer), ecs = Inf),
        ρ = (; cell = fill(1.0, nlayer), ecs = 1.0),
        κ = (;
            cell_interfaces = fill(1e-4, nlayer - 1),
            cell_boundaries = fill(0, nlayer),
            cell_ecs = 1e-4,
            ecs = 0,
        ),
        γ = 2.67513e-4,
    )
end

function get_coeffs(setup::CylinderSetup)
    nlayer = length(setup.groundsetup.layersizes)
    coefficients(
        setup;
        D = (; cell = [0.002 * I(3) for _ = 1:nlayer], ecs = 0.002 * I(3)),
        T₂ = (; cell = fill(Inf, nlayer), ecs = Inf),
        ρ = (; cell = fill(1.0, nlayer), ecs = 1.0),
        κ = (;
            cell_interfaces = fill(1e-4, nlayer - 1),
            cell_boundaries = fill(0, nlayer),
            cell_ecs = 1e-4,
            ecs = 0,
        ),
        γ = 2.67513e-4,
    )
end

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
