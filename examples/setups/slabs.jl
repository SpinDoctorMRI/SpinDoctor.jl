# Floating point type for simulations
T = Float64

groundsetup = PlateSetup{T}(;
    depth = 50.0,
    widths = fill(5.0, 5),
)

# Geometrical setup
setup = ExtrusionSetup(;
    name = "slabs/someslabs",
    groundsetup,
    height = 50.0,
    bend = 0.0,
    twist = 0.0,
    refinement = 10.0,
)

ncell = length(groundsetup.widths)

# Get compartimentalized coefficient vectors
coeffs = coefficients(
    setup;
    D = [0.002 * I(3) for _ = 1:ncell],
    T₂ = fill(Inf, ncell),
    ρ = fill(1.0, ncell),
    κ = (; interfaces = fill(1e-4, ncell - 1), boundaries = fill(0.0, ncell)),
    γ = 2.67513e-4,
)
