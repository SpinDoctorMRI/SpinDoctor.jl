# Geometry and meshes

A simulation starts from a geometry recipe (a "setup"), which [`create_geometry`](@ref)
turns into a finite element mesh in three steps:

1. a configuration of cells (positions and radii), random by default
   (`SpinDoctor.create_cells`);
2. the surfaces of the cells and of the extracellular space (ECS)
   (`SpinDoctor.create_surfaces`);
3. a mesh of triangles (2D, with Triangle) or tetrahedra (3D, with TetGen), split into
   compartments (`SpinDoctor.create_mesh`, `SpinDoctor.split_mesh`).

```julia
setup = CylinderSetup(; ncell = 3, layersizes = [0.6, 1.0], ecs = BoxECS(; margin = 1.0))
mesh, surfaces, cells = create_geometry(setup; rng = Xoshiro(1))
```

See [Geometry recipes](@ref geometries) for pictures of all of them.

## Recipes

| Recipe                  | Dimension | Cells                                         |
| :---------------------- | :-------: | :-------------------------------------------- |
| [`PlateSetup`](@ref)    | 2D        | stacked plates                                |
| [`SlabSetup`](@ref)     | 3D        | stacked slabs (extruded plates)               |
| [`DiskSetup`](@ref)     | 2D        | disks, with layers                            |
| [`CylinderSetup`](@ref) | 3D        | cylinders (extruded disks), with layers       |
| [`SphereSetup`](@ref)   | 3D        | spheres, with layers                          |
| [`VoronoiSetup`](@ref)  | 3D        | polyhedral cells filling a periodic box       |
| [`NeuronSetup`](@ref)   | 3D        | a neuron from a volume mesh (Gmsh format)     |

Disks, cylinders and spheres may have several layers (`layersizes`, relative radii), for
example an axon in a myelin sheath, or a nucleus in a cell. Extruded setups (slabs and
cylinders) may be bent and twisted ([`ExtrusionSetup`](@ref)).

## Extracellular space

The cells of disks, cylinders, spheres and neurons may be immersed in an ECS:

- [`NoECS`](@ref): no ECS;
- [`BoxECS`](@ref): a box around the cells, with a margin;
- [`ConvexHullECS`](@ref): the convex hull of the cells, with a margin;
- [`PeriodicBoxECS`](@ref): a sample of periodic tissue (disks and cylinders).

## Samples of periodic tissue

For dense tissue, [`PeriodicBoxECS`](@ref) and [`VoronoiSetup`](@ref) fill a box with
cells, cut by its faces, with their periodic images entering on the opposite faces:

- [`periodic_packing`](@ref) packs disks or spheres with gaps between them, up to about
  80 % of the area or 60 % of the volume (counting the gaps);
- [`VoronoiSetup`](@ref) creates polyhedral cells for any volume fraction.

The meshes are periodic (the opposite faces of the box have the same points), so that the
solvers can use periodic boundary conditions (see [Solvers](@ref solvers)). Then the
sample represents a larger piece of tissue, without walls.

## Cell configurations

The cell configuration `cells` returned by [`create_geometry`](@ref) may be given back to
create the same cells with other parameters (e.g. a finer mesh), or replaced by a
configuration of your own, for example from another packing algorithm:

```julia
mesh, = create_geometry(setup; cells)
```

The random configurations are drawn from the random number generator `rng`.

## Compartments and boundaries

The mesh is split into compartments, numbered as follows:

- plates and slabs: one per plate;
- disks, cylinders and spheres: the inner layers of all the cells, then the next layers,
  ..., then the ECS (e.g. with two layers and three cells: the axons 1-3, the myelin 4-6,
  the ECS 7);
- neurons: the neuron, then the ECS.

The boundaries are the interfaces between pairs of compartments, followed by the outer
boundaries of each compartment. The coefficients of the recipes ([`coefficients`](@ref))
follow this numbering. [`compartment_mesh`](@ref) extracts one compartment as a mesh of its
own, and [`box_faces`](@ref) finds the facets of a box-shaped compartment on each face.

## Files

By default, [`create_geometry`](@ref) writes no files. With `using JLD2`, the geometry can
be stored in a file and loaded again when it was created for the same setup:

```julia
using JLD2
mesh, surfaces, cells = create_geometry(setup; cache = "geometry.jld2")
```

With `savedir`, the geometry is written as text files in the TetGen format of the MATLAB
version of SpinDoctor.
