# Floating point type for simulations
T = Float64
ncell = 5
# Geometrical setup
setup = PlateSetup{T}(;
    name = "plates/someplates",
    width = 50.0,
    depth = 50.0,
    heights = fill(5.0, ncell),
    deform_angle = (; bend = 0.0, twist =  0.0),
    refinement = 10.0,
    D = [0.002 * I(3) for _ = 1:ncell],
    T₂ = fill(Inf, ncell),
    ρ = fill(1.0, ncell),
    κ = (; interfaces = fill(1e-4, ncell - 1), boundaries = fill(0.0, ncell)),
)