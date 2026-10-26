# Getting started

## Installation

SpinDoctor requires Julia 1.10 or later. In the Julia REPL, press `]` to enter the package
mode and run

```
add SpinDoctor
```

## A first simulation

The signal of three myelinated axons in an extracellular space, for a pulsed gradient
across the axons:

```julia
using SpinDoctor
using LinearAlgebra

# Geometry: three cylinders with two layers (axon and myelin) in an ECS
setup = CylinderSetup(;
    ncell = 3,
    rmin = 2.0,
    rmax = 4.0,
    layersizes = [0.6, 1.0],
    ecs = ConvexHullECS(; margin = 1.0),
    height = 20.0,
)
mesh, surfaces, cells = create_geometry(setup)

# Coefficients of the layers and the ECS (µm, µs)
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

# PGSE sequence (δ = 5 ms, Δ = 20 ms) with b = 1000 s/mm²
gradient = ScalarGradient([1.0, 0.0, 0.0], PGSE(5000.0, 20000.0); bvalue = 1000, gamma = model.gamma)

u = solve_btpde(model, matrices, gradient)    # Magnetization at the echo time
signal = compute_signal(matrices.M, u)        # Signal
adcs = solve_hadc(model, matrices, gradient)  # ADC of each compartment
```

The [Myelinated axons](@ref quickstart) example explains this code
step by step, with figures. The units are explained in [Models and gradients](@ref).

## Optional packages

Some features are loaded with other packages:

| Feature                                      | Packages                                    |
| :------------------------------------------- | :------------------------------------------ |
| Plotting ([Plotting and output](@ref))       | a Makie backend: `CairoMakie`, `GLMakie`, `WGLMakie` |
| Geometry files (`create_geometry(setup; cache)`) | `JLD2`                                  |
| Adaptive time stepping (e.g. `QNDF()`)       | an OrdinaryDiffEq solver package, e.g. `OrdinaryDiffEqBDF` |
| NVIDIA GPUs ([GPU](@ref))                    | `CUDA`, `CUDSS`                             |

The matrix assembly uses all Julia threads: start Julia with `--threads=auto`.

The main functions and types are exported. Lower level building blocks, such as the mesh
type [`SpinDoctor.FEMesh`](@ref) or [`SpinDoctor.assemble_mass_matrix`](@ref), are public
but not exported (see the [API reference](@ref)).

## Examples

The [examples](examples/index.md) are scripts in the `examples` folder of the repository, with
their own environment:

```bash
julia --project=examples examples/quickstart.jl
```
