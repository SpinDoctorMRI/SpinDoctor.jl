# [Solvers](@id solvers)

## Overview

| Method                              | Function                                    | Gives                         | Assumptions                         |
| :---------------------------------- | :------------------------------------------ | :---------------------------- | :---------------------------------- |
| Bloch-Torrey equation               | [`solve_btpde`](@ref)                       | magnetization, signal         | none (reference)                    |
| Homogenized ADC                     | [`solve_hadc`](@ref), [`solve_hadc_tensors`](@ref) | ADC, apparent diffusion tensor | impermeable membranes, low b  |
| Matrix formalism                    | [`solve_laplace`](@ref), [`solve_mf`](@ref) | magnetization, signal         | truncated eigenfunction basis       |
| Matrix formalism, effective tensor  | [`compute_mf_diffusion_tensor`](@ref)       | apparent diffusion tensor     | low b                               |
| Kärger model                        | [`solve_karger`](@ref)                      | signal                        | well mixed compartments             |
| Short diffusion time approximation  | [`compute_adc_sta`](@ref)                   | ADC                           | short diffusion times               |
| Analytical matrix formalism         | [`solve_analytical_mf`](@ref)               | signal                        | layered plates, cylinders, spheres  |

The ADC follows from the signals of several b-values with [`fit_adc`](@ref), and apparent
diffusion tensors from ADCs in several directions with [`fit_tensors`](@ref). See
[Apparent diffusion coefficients](@ref adc-methods) for a comparison.

## Time stepping

[`solve_btpde`](@ref) uses a [`ThetaSolver`](@ref) by default (Crank-Nicolson): one sparse
factorization per interval for PGSE-like sequences, and GMRES with a preconditioner for
other gradients. The time step is a thousandth of the echo time by default:

```julia
u = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 50.0))
signal = compute_signal(matrices.M, u)
```

The adaptive ODE solvers of OrdinaryDiffEq are also available, after loading one of its
solver packages, for example `using OrdinaryDiffEqBDF`:

```julia
u = solve_btpde(model, matrices, gradient, QNDF(); abstol = 1e-8, reltol = 1e-6)
```

## Periodic boundary conditions

The meshes of the samples of periodic tissue ([`PeriodicBoxECS`](@ref),
[`VoronoiSetup`](@ref)) are periodic. With

```julia
matrices = assemble_matrices(model; periodic = true)
```

the points on opposite faces of the box are identified, and [`solve_btpde`](@ref) (with
[`ThetaSolver`](@ref)) and [`solve_hadc`](@ref) use periodic boundary conditions instead of
impermeable box faces (see [Periodic boundary conditions](@ref periodic-theory)). The
signals and ADCs then do not depend on the position of the cells in the box, and the water
of the ECS is not trapped by the box. See [Dense tissue with periodic boundary conditions](@ref periodic) and
[White matter](@ref white-matter).

## Callbacks

[`solve_btpde`](@ref) calls the `callbacks` after each time step: [`Printer`](@ref) prints
information, [`VTKWriter`](@ref) saves the magnetization for
[ParaView](https://www.paraview.org), and [`Plotter`](@ref) plots it (see
[Plotting and output](@ref)).
