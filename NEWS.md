# Changes

## v0.2.0

This release runs on current Julia versions (1.10 and later) with current versions of all
dependencies. The core no longer depends on the SciML stack, Arpack or Expokit. It includes the modernization of
[#47](https://github.com/SpinDoctorMRI/SpinDoctor.jl/pull/47) by Morgan G. Hough (Makie
extension, split OrdinaryDiffEq packages).

### Breaking changes

- The problem types and `solve` methods are replaced by solver functions:
  `solve_btpde(model, matrices, gradient[, solver])`, `solve_hadc(model, matrices,
  gradient)`, `solve_karger(model, difftensors, gradient; timestep)`,
  `solve_laplace(model, matrices; neig_max)`, `solve_mf(model, matrices, lap_eig,
  gradient)`, `solve_analytical_laplace(; coeffs...)` and `solve_analytical_mf(laplace,
  volumes, gradient)`.
- Plotting requires loading a Makie backend (package extension).
- Time stepping is built in: `ThetaSolver(; theta, timestep)` for the BTPDE (any gradient;
  one factorization per interval for interval-wise constant gradients, preconditioned GMRES
  otherwise) replaces `IntervalConstantSolver`, and the HADC uses Crank-Nicolson.
  OrdinaryDiffEq is an optional extension (`using OrdinaryDiffEqBDF; solve_btpde(model,
  matrices, gradient, QNDF())`), and `QNDF` is no longer exported. `solve_hadc` takes a
  `timestep` instead of an ODE solver and tolerances.
- `create_geometry` writes no files by default, and reuses existing files only if they were
  created for the same setup (`recreate = false` by default). Geometries can be cached in
  JLD2 files (`cache`, with `using JLD2`).
- `NeuronSetup` takes the Gmsh volume mesh as `meshfile`.
- `assemble_mass_matrix` takes `elementweights` or `nodeweights` as keywords.
- `TightWrapECS`, which was not implemented, is removed.
- ASCII names in the API: the coefficients and model fields are `rho`, `D`, `T2`, `kappa`
  and `gamma`, the time profiles have the fields `delta` and `Delta`, `int_F²` is
  `int_F2`, and the analytical permeabilities `W` are `kappa`.
- Renamed: `get_cmpt_volumes` → `compartment_volumes`, `get_mesh_volumes` →
  `element_geometry`, `get_mesh_surface` and `get_mesh_surfacenormals` → `facet_geometry`
  (outward normals when the elements are given), `couple_flux_matrix` →
  `assemble_flux_matrix`, `savefield(mesh, u, filename)` → `save_field(filename, mesh, u)`,
  `create_fibonacci_sphere` → `fibonacci_sphere`. The text geometry writers take the file
  name first (`save_cells`, `save_surfaces`, `save_mesh`).
- Only the main API is exported. Building blocks (e.g. `FEMesh`, `assemble_mass_matrix`,
  `integral`, `int_F2`, `intervals`, `echotime`) are public but not exported.

### Fixed

- Outward normals of non-convex compartments (ECS, neurons, bent cells) were partly
  inverted, which made the HADC of these compartments wrong.
- The mass matrix of line segments was twice too large, which doubled the permeability
  coupling (and the HADC surface terms) in 2D.
- Boundary facets with all their points on an interface were also assigned to the
  neighboring compartment, so that layered or bent cylinders could not be assembled.
- The numerical `integral` and `int_F²` of custom time profiles returned tuples.
- The Laplace eigenproblem did not converge for many eigenvalues; it now uses
  shift-and-invert Lanczos iterations (KrylovKit instead of Arpack).
- The analytical matrix formalism missed eigenvalues (closely spaced ones for weakly
  coupled layers, and high angular orders), so that its signal did not converge with the
  number of modes (e.g. for four permeable plates: 0.66, 0.76, 0.65, 0.71 for decreasing
  length scales; now 0.6526 for all). Bessel overflows for high orders are handled.
- 2D VTK output, ignored keywords of `unitsphere`, the precision of the text geometry files,
  the quadratic-time extraction of neuron boundaries, `Printer` ignoring `nupdate`, a
  division by zero in the theta-scheme for intervals shorter than half a time step,
  and `ScalarGradient` and `GeneralGradient` not being `AbstractGradient`s.

### Added

- Multithreaded matrix assembly with static element kernels, much faster and with far
  less memory than before.
- `compartment_mesh`, `box_faces`, `compartment_weights`, given cell configurations
  (`create_geometry(setup; cells)`), random number generators (`rng`), and
  `save_geometry`/`load_geometry`.
- `solve_hadc_tensors`: the apparent diffusion tensor of each compartment from the HADC
  model, for all directions at once (one solve per axis).
- `ScalarGradient(dir, profile; bvalue, gamma)` and `bvalue(gradient, gamma)`; time
  profiles accept integer arguments.
- `fit_adc(bvalues, signals; degree)`.
- Dense samples of periodic tissue: `periodic_packing` packs polydisperse disks or spheres
  with gaps in a periodic box (collective rearrangement), `PeriodicBoxECS` fills the box
  with disks or cylinders cut by its faces (without the free ECS layer of `BoxECS`), and
  `VoronoiSetup` creates polyhedral cells (shrunk Voronoi cells) for any cell volume
  fraction. The TetGen meshing accepts polygon facets and several seeds per region.
- Time stepping on NVIDIA GPUs with the sparse direct solvers of cuDSS (`using CUDA, CUDSS`):
  `ThetaSolver(; arraytype = CuArray)` for the BTPDE, and `solve_hadc(...; arraytype =
  CuArray)`. For a mesh of 7e4 nodes, an RTX 4090 is 10 to 18 times faster than 8 CPU
  threads. The GPU tests have their own environment in `test/gpu`.
- The HADC solves all compartments together, and the BTPDE with non-constant gradients
  applies the matrix of the gradient without forming it at every step.
- Value-based tests (analytical solutions, exact identities), and regression outputs of
  v0.1.0 in `test/regression`.
