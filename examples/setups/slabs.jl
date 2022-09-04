# Name for saving meshfiles and data
name = "slabs"
meshdir = joinpath("meshfiles", name)

# Geometrical setup
ncell = 5
setup = SlabSetup(;
    depth = 50.0,
    widths = fill(5.0, ncell),
    height = 50.0,
    bend = 0.03,
    twist = π / 6,
    refinement = 10.0,
)

# Get compartimentalized coefficient vectors
coeffs = coefficients(
    setup;
    D = [0.002 * I(3) for _ = 1:ncell],
    T₂ = fill(Inf, ncell),
    ρ = fill(1.0, ncell),
    κ = (; interfaces = fill(1e-4, ncell - 1), boundaries = fill(0.0, ncell)),
    γ = 2.67513e-4,
)
