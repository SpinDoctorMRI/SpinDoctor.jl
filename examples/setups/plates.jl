# Floating point type for simulations
T = Float64

# Name for saving meshfiles and data
name = "plates"
savedir = joinpath("data", name)

# Geometrical setup
setup = PlateSetup(; depth = 20.0, widths = fill(5.0, 10), refinement = 1.0)

ncell = length(setup.widths)

# Get compartimentalized coefficient vectors
coeffs = coefficients(
    setup;
    D = [0.002 * I(2) for _ = 1:ncell],
    T2 = fill(Inf, ncell),
    rho = fill(1.0, ncell),
    kappa = (; interfaces = fill(1e-4, ncell - 1), boundaries = fill(0.0, ncell)),
    gamma = 2.67513e-4,
)
