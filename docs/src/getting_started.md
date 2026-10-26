# Getting started

## Installation

To install SpinDoctor, run the following commands:

```
using Pkg
Pkg.add("SpinDoctor")
```

Multiple examples are found in the `examples` directory, with their own environment
(`julia --project=examples`).

Some features require loading additional packages:

- plotting ([`plot_mesh`](@ref), [`plot_field`](@ref), the [`Plotter`](@ref) callback)
  requires a Makie backend, e.g. `using CairoMakie` or `using GLMakie`;
- saving geometries ([`save_geometry`](@ref), `create_geometry(setup; cache)`) requires
  `using JLD2`.

The finite element matrices are assembled in parallel on all Julia threads (start Julia
with `--threads=auto`).
