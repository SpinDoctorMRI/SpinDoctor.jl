````@raw html
---
layout: home

hero:
  name: SpinDoctor
  text: Diffusion MRI simulations
  tagline: The Bloch-Torrey equation on finite element meshes of cells and tissue, in Julia
  image:
    src: /logo.svg
    alt: SpinDoctor
  actions:
    - theme: brand
      text: Getting started
      link: /getting_started
    - theme: alt
      text: Examples
      link: /examples
    - theme: alt
      text: View on GitHub
      link: https://github.com/SpinDoctorMRI/SpinDoctor.jl

features:
  - icon: 🧠
    title: Cells and tissue
    details: Recipes for plates, disks, cylinders, spheres and dense samples of periodic tissue, with layers, an extracellular space and permeable membranes, and realistic neurons.
    link: /examples/generated/geometries
  - icon: 🧲
    title: Signals and ADCs
    details: Solve the Bloch-Torrey equation for any gradient sequence, or compute the ADC with the homogenized model, the matrix formalism, the Kärger model and short-time approximations.
    link: /examples/generated/quickstart
  - icon: 🔁
    title: Periodic boundaries
    details: Samples of densely packed tissue with periodic boundary conditions, free of the walls of a finite sample.
    link: /examples/generated/white_matter
  - icon: 🚀
    title: Fast
    details: Multithreaded assembly, built-in time stepping with sparse direct solvers, and NVIDIA GPUs through CUDA and cuDSS.
    link: /manual/gpu
---
````

## What is SpinDoctor?

SpinDoctor simulates diffusion magnetic resonance imaging (dMRI) in biological tissue
[Li2019](@cite). The water in the tissue diffuses in compartments (cells, their layers,
and the extracellular space) separated by permeable membranes, and the magnetization of its
spins follows the Bloch-Torrey partial differential equation. SpinDoctor solves this
equation with P1 finite elements on meshes created by
[TetGen](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) [Si2015](@cite) and
[Triangle](https://www.cs.cmu.edu/~quake/triangle.html), and computes:

- the dMRI signal for scalar (PGSE, double PGSE, OGSE) and arbitrary three-dimensional
  gradient sequences;
- apparent diffusion coefficients (ADC): by fitting the signal, with the homogenized ADC
  model, with the matrix formalism (Laplace eigenfunctions), and with the short diffusion
  time approximation;
- analytical solutions for layered plates, cylinders and spheres, for validation.

The compartments may have different initial spin densities, diffusion tensors, T2
relaxation and permeabilities [Nguyen2014d](@cite).

This is the Julia version of the [SpinDoctor](https://github.com/SpinDoctorMRI/SpinDoctor)
MATLAB toolbox.

## Citing

If you use SpinDoctor in your research, please cite [Li2019](@cite), and
[Agdestein2021](@cite) for the matrix formalism.
