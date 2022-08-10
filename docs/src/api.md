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
Model
initial_conditions
```

## Postprocessing

```@docs
savefield
compute_adc_sta
fit_adc
fit_tensors
plot_mesh
plot_field
```

## Problems

SpinDoctor allows for considering a wide range of diffusion MRI problems. These can be
solved for using the [`solve`](@ref) function.

```@docs
BTPDE
HADC
Karger
Laplace
MatrixFormalism
AnalyticalLaplace
AnalyticalMatrixFormalism
IntervalConstantSolver
solve
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
SphereSetup
NeuronSetup
coefficients
analytical_coefficients
create_geometry
```
