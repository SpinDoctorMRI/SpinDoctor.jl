# SpinDoctor

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://agdestein.github.io/SpinDoctor.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://agdestein.github.io/SpinDoctor.jl/dev)
[![Build Status](https://github.com/agdestein/SpinDoctor.jl/workflows/CI/badge.svg)](https://github.com/agdestein/SpinDoctor.jl/actions)
[![Coverage](https://codecov.io/gh/agdestein/SpinDoctor.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/agdestein/SpinDoctor.jl)

This is a Julia implementation of the SpinDoctor toolbox. The original MATLAB toolbox can be found at https://github.com/jingrebeccali/SpinDoctor.

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.

SpinDoctor can be used

1) to solve the Bloch-Torrey partial differential equation (BTDPE) to obtain the dMRI signal (the toolbox provides a way of robustly fitting the dMRI signal to obtain the fitted Apparent Diffusion Coefficient (ADC));
2) to solve the diffusion equation for the homogenized ADC (HADC) model to obtain the ADC;
3) a short-time approximation formula for the ADC is also included in the toolbox for comparison with the simulated ADC;
4) Compute the dMRI signal using a matrix formalism (MF) analytical solution based Laplace eigenfunctions.

The PDEs and Laplace eigenvalue decompositions are solved by P1 finite elements combined with built-in MATLAB routines for solving ordinary differential equations.
The finite element mesh generation is performed using an external package called TetGen.

SpinDoctor has support for the following features:
1. multiple compartments with compartment-wise constant
	* initial spin densities,
	* diffusion coefficients or diffusion tensors, and
	* T2-relaxation coefficients;
2. permeable membranes between compartments for the BTPDE and MF (the HADC assumes negligible permeability);
3. built-in diffusion-encoding pulse sequences, including
	* the Pulsed Gradient Spin Echo (PGSE) and double-PGSE,
	* the Ocsillating Gradient Spin Echo (OGSE), or
	* custom pulse sequences;
4. uniformly distributed gradient directions in 2D and 3D for high angular resolution diffusion imaging (HARDI)

SpinDoctor also comes with a geometry generation module, allowing for
1. spherical cells with a nucleus;
2. cylindrical cells with a myelin layer;
3. an extra-cellular space (ECS) enclosed in either
	* a box,
	* a convex hull, or
	* a tight wrapping around the cells;
4. deformation of canonical cells by bending and twisting.

In addition, a variety of neuron meshes is available, whose surface geometries were extracted from [NeuroMopho.org](http://neuromorpho.org). The neurons may also be enclosed in an extracellular space as described above.


### Spinning spindle spins in SpinDoctor

![Spindle](misc/spindle.gif)

The above graphic visualizes the magnetization as a z-displacement for the spindle neuron geometry `03b_spindle4aACC` (extracted from NeuroMorpho). The gradient is a PGSE sequence in the x-direction. The magnetization was computed for 200 time steps, and the exported vtk sequence was visualized using [Paraview](https://www.paraview.org).

 
## Getting started

1) The base folder contains a commented general purpose driver `driver.jl`.
2) The input files for the drivers are found in the folder `setups`, and define the structures needed for the simulations.
3) Multiple neuron meshes are found in the folder `mesh_files`. These can be loaded in the file `setups/neuron.jl`.
4) The user guide is found [here](https://github.com/jingrebeccali/SpinDoctor/blob/master/user_guide.pdf).


## How to cite us

The paper about SpinDoctor can be found at https://arxiv.org/abs/1902.01025.

If you use our software for research, please consider citing us:

```bibtex
@article{Li2019,
  author  = {Jing-Rebecca Li and Van-Dang Nguyen and Try Nguyen Tran and Jan Valdman and Cong-Bang Trang and Khieu Van Nguyen and Duc Thach Son Vu and Hoang An Tran and Hoang Trong An Tran and Thi Minh Phuong Nguyen},
  doi     = {https://doi.org/10.1016/j.neuroimage.2019.116120},
  issn    = {1053-8119},
  journal = {NeuroImage},
  pages   = {116120},
  title   = {{SpinDoctor: A MATLAB toolbox for diffusion MRI simulation}},
  url     = {http://www.sciencedirect.com/science/article/pii/S1053811919307116},
  volume  = {202},
  year    = {2019}
}
```
