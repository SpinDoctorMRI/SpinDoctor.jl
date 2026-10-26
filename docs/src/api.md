# API Reference

## Gradients

These functions deal with gradient pulse sequences.

```@docs
GeneralGradient
ScalarGradient
PGSE
DoublePGSE
CosOGSE
SinOGSE
integral
int_F²
intervals
isconstant
echotime
```

## Geometry

```@docs
assemble_matrices
split_mesh
split_field
get_cmpt_volumes
compartment_mesh
box_faces
Model
initial_conditions
```

## Postprocessing

```@docs
savefield
compute_adc_sta
fit_adc
fit_tensors
compartment_weights
```

## Plotting

The plotting functions require a Makie backend, e.g. `using CairoMakie` or
`using GLMakie`.

```@docs
plot_surfaces
plot_mesh
plot_field
plot_hardi
```

## Solvers

SpinDoctor allows for considering a wide range of diffusion MRI problems, each with its own
solver function.

```@docs
solve_btpde
IntervalConstantSolver
solve_hadc
solve_karger
solve_laplace
solve_mf
solve_analytical_laplace
solve_analytical_mf
```

## Matrix formalism

```@docs
length2eig
eig2length
limit_lengthscale
compute_mf_diffusion_tensor
```

## Utils

```@docs
unitcircle
unitsphere
compute_signal
```

## Callbacks

Callbacks are called after every time step when solving the BTPDE.

```@docs
Printer
VTKWriter
Plotter
```

## Recipes

SpinDoctor comes with pre-configures recipes for creating finite element meshes of plates,
cylinders, spheres, and neurons.

```@docs
PlateSetup
DiskSetup
ExtrusionSetup
CylinderSetup
SlabSetup
SphereSetup
NeuronSetup
NoECS
BoxECS
ConvexHullECS
coefficients
analytical_coefficients
create_geometry
```
