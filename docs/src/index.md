```@meta
CurrentModule = SpinDoctor
```

# SpinDoctor

This is the documentation for [SpinDoctor](https://github.com/agdestein/SpinDoctor.jl), and
is based on
[the SpinDoctor User Guide](https://github.com/jingrebeccali/SpinDoctor/blob/master/user_guide.pdf).

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic
resonance imaging (dMRI) for prototyping purposes.

SpinDoctor can be used

1. to solve the Bloch-Torrey partial differential equation (BTDPE) to obtain the dMRI signal
   (the toolbox provides a way of robustly fitting the dMRI signal to obtain the fitted
   Apparent Diffusion Coefficient (ADC));
2. to solve the diffusion equation for the homogenized ADC (HADC) model to obtain the ADC;
3. a short-time approximation formula for the ADC is also included in the toolbox for
   comparison with the simulated ADC;
4. Compute the dMRI signal using a matrix formalism (MF) analytical solution based Laplace
   eigenfunctions.

The PDEs and Laplace eigenvalue decompositions are solved by P1 finite elements. The
geometry recipes create surface triangulations that are passed to
[TetGen](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) to perform the finite
element mesh generation.

SpinDoctor has support for the following features:
1. multiple compartments connected through permeable membranes, with different
    * initial spin densities,
    * diffusion tensors,
    * T2-relaxation, and
    * permeability coefficients for the BTPDE and MF (the HADC assumes negligible
      permeability);
3. diffusion-encoding gradient pulse sequences, including
    * the pulsed gradient spin echo sequence (PGSE),
    * the double-PGSE,
    * the oscillating gradient spin echo (OGSE), and
    * custom three-dimensional pulse sequences `g(t) = (gx(t), gy(t), gz(t))`;
4. uniformly distributed gradient directions in 2D and 3D for high angular resolution
   diffusion imaging (HARDI).

SpinDoctor also comes with a geometry generation module, allowing for

1. spherical cells with a nucleus;
2. cylindrical cells with a myelin layer;
3. an extra-cellular space (ECS) enclosed in either a box, a convex hull, or a tight
   wrapping around the cells;
4. deformation of canonical cells by bending and twisting.

In addition, a variety of neuron meshes is available, whose surface geometries were
extracted from [NeuroMopho.org](http://neuromorpho.org). The neurons may also be enclosed in
an extracellular space as described above.
