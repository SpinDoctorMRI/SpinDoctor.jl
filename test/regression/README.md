# Regression outputs

Numerical outputs of SpinDoctor on a set of small setups, used to check that
refactoring does not change results unintentionally.

- `record.jl` computes the outputs (matrices, ADCs, signals, eigenvalues) for
  each setup and writes them to a TOML file.
- `meshes/` holds the geometries, read back by `record.jl` so that every
  version starts from identical meshes.
- `baseline.toml` holds the outputs of SpinDoctor v0.1.0 (commit `f824161`,
  January 2023), run on Julia 1.8.5 with the environment in `baseline-env/`.
  It was recorded with the version of `record.jl` from commit `7cf8805`,
  which uses the API of that time (`solve(BTPDE(; model, matrices), ...)`).
  It is never changed.
- `reference.toml` holds the currently accepted outputs. It starts as a copy
  of the baseline. When a bug fix changes outputs on purpose, it is
  regenerated and the change is listed below.
- `compare.jl` compares two output files with per-quantity tolerances.

To check the current code:

```bash
julia --project=test/regression test/regression/record.jl /tmp/outputs.toml
julia test/regression/compare.jl /tmp/outputs.toml test/regression/reference.toml
```

## Failures in the baseline

The baseline records an error instead of values for:

- `cylinders_bent_hull` and `cylinders_layers_hull`: `assemble_matrices`
  fails, since `split_mesh` assigns boundaries to the wrong compartments.
- `neuron_box` and `plates_2d`: the HADC fails inside the linear solver of
  the OrdinaryDiffEq version of that time.

## Intended changes to the reference

1. Segment mass matrix (commit "Fix the P1 mass matrix of line segments"): in
   2D (`plates_2d`, `disks_box_2d`), `matrix_Q_*` and `G_*` are halved, and
   all signals, eigenvalues and ADCs change accordingly (by 5–10 % for the
   signals here).
2. Outward normals (commit "Fix the outward normals of non-convex
   compartments"): `G_*`, `hadc` and `karger` change for the ECS of all
   setups with an ECS, for the neuron, and for the bent cylinders. The other
   compartments are unchanged.
3. Bent and layered cylinders (commit "Assign boundary facets to
   compartments by element faces"): these setups now run; their outputs are
   recorded for the first time.

4. Analytical eigenvalues (commit "Hand-written time stepping..."): `analytical`
   changes by 1.5 %, since eigenvalues were missed.
5. Time stepping (same commit): `hadc`, `hadc_ogse` and
   `karger` are computed with the Crank-Nicolson HADC instead of QNDF, and change within
   the tolerance of QNDF (at most 2e-4). The BTPDE keys `btpde_qndf` and `btpde_general`
   still use QNDF, through the OrdinaryDiffEq extension.

The current `reference.toml` was recorded after these changes, on Julia 1.12, with the
environment `Project.toml` of this folder.
