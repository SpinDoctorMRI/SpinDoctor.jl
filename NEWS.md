# Changes

## v0.2.0

This release runs on current Julia versions (1.10 and later) with current versions of all
dependencies. It includes the modernization of
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
- `create_geometry` writes no files by default, and reuses existing files only if they were
  created for the same setup (`recreate = false` by default). Geometries can be cached in
  JLD2 files (`cache`, with `using JLD2`).
- `NeuronSetup` takes the Gmsh volume mesh as `meshfile`.
- `assemble_mass_matrix` takes `elementweights` or `nodeweights` as keywords.
- `TightWrapECS`, which was not implemented, is removed.
- ASCII names in the API: the coefficients and model fields are `rho`, `D`, `T2`, `kappa`
  and `gamma`, the time profiles have the fields `delta` and `Delta`, `int_F²` is
  `int_F2`, the analytical permeabilities `W` are `kappa`, and `IntervalConstantSolver`
  takes `theta`.
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
- 2D VTK output, ignored keywords of `unitsphere`, the precision of the text geometry files,
  the quadratic-time extraction of neuron boundaries, `Printer` ignoring `nupdate`, a
  division by zero in `IntervalConstantSolver` for intervals shorter than half a time step,
  and `ScalarGradient` and `GeneralGradient` not being `AbstractGradient`s.

### Added

- Multithreaded matrix assembly with static element kernels, much faster and with far
  less memory than before.
- `compartment_mesh`, `box_faces`, `compartment_weights`, given cell configurations
  (`create_geometry(setup; cells)`), random number generators (`rng`), and
  `save_geometry`/`load_geometry`.
- `ScalarGradient(dir, profile; bvalue, gamma)` and `bvalue(gradient, gamma)`; time
  profiles accept integer arguments.
- `fit_adc(bvalues, signals; degree)`.
- Value-based tests (analytical solutions, exact identities), and regression outputs of
  v0.1.0 in `test/regression`.
