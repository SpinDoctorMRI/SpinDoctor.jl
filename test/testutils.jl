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
        T2 = fill(Inf, nlayer),
        rho = fill(1.0, nlayer),
        kappa = (; interfaces = fill(1e-4, nlayer), boundaries = fill(0.0, nlayer)),
        gamma = 2.67513e-4,
    )
end

function get_coeffs(setup::SlabSetup)
    nlayer = length(setup.groundsetup.widths)
    coefficients(
        setup;
        D = [0.002 * I(3) for _ = 1:nlayer],
        T2 = fill(Inf, nlayer),
        rho = fill(1.0, nlayer),
        kappa = (; interfaces = fill(1e-4, nlayer), boundaries = fill(0.0, nlayer)),
        gamma = 2.67513e-4,
    )
end

function get_coeffs(setup::DiskSetup)
    nlayer = length(setup.layersizes)
    coefficients(
        setup;
        D = (; cell = [0.002 * I(2) for _ = 1:nlayer], ecs = 0.002 * I(2)),
        T2 = (; cell = fill(Inf, nlayer), ecs = Inf),
        rho = (; cell = fill(1.0, nlayer), ecs = 1.0),
        kappa = (;
            cell_interfaces = fill(1e-4, nlayer - 1),
            cell_boundaries = fill(0, nlayer),
            cell_ecs = 1e-4,
            ecs = 0,
        ),
        gamma = 2.67513e-4,
    )
end

function get_coeffs(setup::CylinderSetup)
    nlayer = length(setup.groundsetup.layersizes)
    coefficients(
        setup;
        D = (; cell = [0.002 * I(3) for _ = 1:nlayer], ecs = 0.002 * I(3)),
        T2 = (; cell = fill(Inf, nlayer), ecs = Inf),
        rho = (; cell = fill(1.0, nlayer), ecs = 1.0),
        kappa = (;
            cell_interfaces = fill(1e-4, nlayer - 1),
            cell_boundaries = fill(0, nlayer),
            cell_ecs = 1e-4,
            ecs = 0,
        ),
        gamma = 2.67513e-4,
    )
end

function get_coeffs(setup::SphereSetup)
    nlayer = length(setup.layersizes)
    coefficients(
        setup;
        D = (; cell = [0.002 * I(3) for _ = 1:nlayer], ecs = 0.002 * I(3)),
        T2 = (; cell = fill(Inf, nlayer), ecs = Inf),
        rho = (; cell = fill(1.0, nlayer), ecs = 1.0),
        kappa = (;
            cell_interfaces = fill(1e-4, nlayer - 1),
            cell_boundaries = fill(0, nlayer),
            cell_ecs = 1e-4,
            ecs = 0,
        ),
        gamma = 2.67513e-4,
    )
end

get_coeffs(setup::NeuronSetup) = coefficients(
    setup;
    D = (; neuron = 0.002 * I(3), ecs = 0.002 * I(3)),
    T2 = (; neuron = Inf, ecs = Inf),
    rho = (; n = 1.0, neuron = 1.0, ecs = 1.0),
    kappa = (; neuron_ecs = 1e-4, neuron = 0.0, ecs = 0.0),
    gamma = 2.67513e-4,
)

get_coeffs(setup::VoronoiSetup) =
    get_coeffs(SphereSetup(; ncell = setup.ncell, ecs = BoxECS()))
