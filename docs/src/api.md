# API Reference

```@docs
SpinDoctor
```

The main functions and types are exported. The building blocks listed under
[Public, not exported](@ref) are used with the module name, e.g. `SpinDoctor.FEMesh`.

## Setup recipes and geometry

SpinDoctor comes with recipes for creating finite element meshes of plates, slabs, disks,
cylinders, spheres and neurons, possibly immersed in an extracellular space (ECS).

```@docs
PlateSetup
SlabSetup
DiskSetup
CylinderSetup
SphereSetup
ExtrusionSetup
NeuronSetup
NoECS
BoxECS
ConvexHullECS
create_geometry
save_geometry
load_geometry
compartment_mesh
box_faces
compartment_volumes
```

## Model and matrices

```@docs
coefficients
analytical_coefficients
Model
initial_conditions
assemble_matrices
```

## Gradients

```@docs
PGSE
DoublePGSE
CosOGSE
SinOGSE
ScalarGradient
GeneralGradient
bvalue
unitcircle
unitsphere
```

## Solvers

SpinDoctor allows for considering a wide range of diffusion MRI problems, each with its own
solver function.

```@docs
solve_btpde
ThetaSolver
solve_hadc
solve_hadc_tensors
compute_adc_sta
solve_karger
solve_laplace
solve_mf
eig2length
length2eig
limit_lengthscale
compute_mf_diffusion_tensor
solve_analytical_laplace
solve_analytical_mf
```

## Post-processing and output

```@docs
compute_signal
split_field
compartment_weights
fit_adc
fit_tensors
save_field
```

## Callbacks

Callbacks are called after every time step when solving the BTPDE.

```@docs
Printer
VTKWriter
Plotter
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

## Public, not exported

### Meshes and finite elements

```@docs
SpinDoctor.FEMesh
SpinDoctor.split_mesh
SpinDoctor.facet_elements
SpinDoctor.element_geometry
SpinDoctor.facet_geometry
SpinDoctor.assemble_mass_matrix
SpinDoctor.assemble_stiffness_matrix
SpinDoctor.assemble_flux_matrix
```

### Geometry generation and files

```@docs
SpinDoctor.AbstractSetup
SpinDoctor.AbstractECS
SpinDoctor.create_cells
SpinDoctor.create_surfaces
SpinDoctor.create_mesh
SpinDoctor.deform_domain!
SpinDoctor.radial_dimension
SpinDoctor.fibonacci_sphere
SpinDoctor.read_gmsh
SpinDoctor.read_cells
SpinDoctor.save_cells
SpinDoctor.read_surfaces
SpinDoctor.save_surfaces
SpinDoctor.read_mesh
SpinDoctor.save_mesh
```

### Gradient sequences

```@docs
SpinDoctor.AbstractTimeProfile
SpinDoctor.AbstractGradient
SpinDoctor.integral
SpinDoctor.int_F2
SpinDoctor.intervals
SpinDoctor.echotime
SpinDoctor.isconstant
SpinDoctor.diffusion_time
SpinDoctor.diffusion_time_sta
SpinDoctor.j_integral
```

### Callbacks

```@docs
SpinDoctor.AbstractCallback
```
