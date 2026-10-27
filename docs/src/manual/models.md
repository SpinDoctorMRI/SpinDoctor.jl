# Models and gradients

## Units

SpinDoctor uses micrometers and microseconds:

| Quantity                     | Unit                | Example                             |
| :--------------------------- | :------------------ | :---------------------------------- |
| Length                       | µm                  |                                     |
| Time                         | µs                  | `PGSE(5000.0, 20000.0)`: 5 and 20 ms |
| Diffusion coefficient        | µm²/µs              | `0.002` is 2 × 10⁻³ mm²/s (water)    |
| b-value                      | µs/µm² (= s/mm²)    | `1000`                              |
| Permeability, relaxivity     | µm/µs (= m/s)       | `1e-5`                              |
| Gradient amplitude           | T/m                 |                                     |
| Gyromagnetic ratio           | rad / (µs µm T/m)   | `2.67513e-4` for water protons      |

## Model

A [`Model`](@ref) contains the mesh and the coefficients of each compartment (initial spin
density `rho`, diffusion tensor `D`, relaxation time `T2`) and of each boundary
(permeability or relaxivity `kappa`), and the gyromagnetic ratio `gamma`.
[`coefficients`](@ref) distributes the coefficients of the layers and the ECS of a recipe
over the compartments and boundaries:

```julia
coeffs = coefficients(
    setup;
    D = (; cell = [0.002 * I(3), 0.0005 * I(3)], ecs = 0.002 * I(3)),
    T2 = (; cell = [Inf, Inf], ecs = Inf),
    rho = (; cell = [1.0, 1.0], ecs = 1.0),
    kappa = (; cell_interfaces = [1e-5], cell_boundaries = [0, 0], cell_ecs = 1e-5, ecs = 0),
    gamma = 2.67513e-4,
)
model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)
```

The diffusion tensor of a compartment may also vary: a vector with one tensor per element
of the compartment (e.g. from diffusion tensor imaging) replaces the constant tensor, in
the same order as `mesh.elements[i]`.

[`assemble_matrices`](@ref) assembles the finite element matrices, in parallel on all Julia
threads (start Julia with `--threads=auto`). [`compartment_weights`](@ref) gives the
contribution of each compartment to the signal without gradient.

## Gradient sequences

A [`ScalarGradient`](@ref) has a fixed direction and a time profile:

- [`PGSE`](@ref): two pulses of duration `delta`, separated by `Delta`;
- [`DoublePGSE`](@ref): two PGSE sequences;
- [`CosOGSE`](@ref) and [`SinOGSE`](@ref): oscillating gradients.

The amplitude is given directly, or by the b-value:

```julia
gradient = ScalarGradient([1.0, 0.0, 0.0], PGSE(5000.0, 20000.0); bvalue = 1000, gamma)
bvalue(gradient, gamma)
```

A [`GeneralGradient`](@ref) is any function of time, with a varying direction. Directions
for high angular resolution imaging are given by [`unitcircle`](@ref) and
[`unitsphere`](@ref).
