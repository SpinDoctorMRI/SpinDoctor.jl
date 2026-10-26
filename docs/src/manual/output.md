# Plotting and output

## Plotting

The plotting functions require a [Makie](https://docs.makie.org) backend: `CairoMakie` for
images, `GLMakie` for interactive windows, or `WGLMakie` for the browser.

- [`plot_mesh`](@ref): the compartments of a mesh, with colors and transparency per
  compartment;
- [`plot_field`](@ref): a field on the mesh, such as the magnetization (magnitude, real or
  imaginary part, or phase);
- [`plot_hardi`](@ref): values (signals, ADCs) in many directions, as a glyph;
- [`plot_surfaces`](@ref): the surfaces of a geometry, before meshing.

The functions with a `!` (`plot_mesh!`, `plot_field!`, `plot_hardi!`) plot into a given
axis, for figures with several panels:

```julia
using CairoMakie
fig = Figure()
plot_mesh!(Axis3(fig[1, 1]; aspect = :data), mesh)
plot_field!(Axis3(fig[1, 2]; aspect = :data), mesh, u)
fig
```

In 3D, the surfaces of the compartments are shown: samples of periodic tissue show the
cross-sections of the cells on the faces of the box.

## Following the time stepping

The [`Plotter`](@ref) callback plots the gradient, the signal attenuation and the
magnetization during the time stepping (live with GLMakie). [`Printer`](@ref) prints
information, and [`VTKWriter`](@ref) saves the magnetization at the time steps for
[ParaView](https://www.paraview.org).

## Files

- [`save_field`](@ref) saves a field on the mesh in the VTK format;
- [`save_geometry`](@ref) and [`load_geometry`](@ref) store a geometry in a JLD2 file
  (`using JLD2`);
- `SpinDoctor.read_gmsh` reads volume meshes in the Gmsh format.
