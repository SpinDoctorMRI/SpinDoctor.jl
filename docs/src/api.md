# API Reference

## Gradients

```@docs
PGSE
CosOGSE
SinOGSE
DoublePGSE
GeneralGradient
ScalarGradient
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

```@docs
IntervalConstanBTPDE
GeneralBTPDE
HADC
Karger
Laplace
MatrixFormalism
AnalyticalLaplace
AnalyticalMatrixFormalism
solve
solve_multigrad
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

```@docs
Printer
VTKWriter
Plotter
```

## Recipes

```@docs
AbstractSetup
PlateSetup
CylinderSetup
SphereSetup
NeuronSetup
coefficients
analytical_coefficients
create_geometry
```
