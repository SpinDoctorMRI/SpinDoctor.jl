# # [Geometry recipes](@id geometries)
#
# SpinDoctor creates finite element meshes from geometry recipes ("setups"): a random
# configuration of cells, the surfaces of the cells and of the extracellular space (ECS),
# and a mesh of triangles (2D) or tetrahedra (3D). This example shows the available
# recipes.

using SpinDoctor
using Random
using CairoMakie

# Colors for the plots: the layers of the cells and the ECS.

axon, myelin, ecs = :steelblue, :orange, (:gray, 0.9)
nothing #hide

# ## Plates and slabs
#
# Stacked plates in 2D ([`PlateSetup`](@ref)), with the given widths, and their extrusion
# to slabs in 3D ([`SlabSetup`](@ref)), here twisted around the vertical axis.

plates = PlateSetup(; depth = 20.0, widths = [3.0, 5.0, 4.0, 6.0], refinement = 2.0)
mesh_plates, = create_geometry(plates)
plot_mesh(mesh_plates)

#-

slabs = SlabSetup(; depth = 10.0, widths = fill(2.5, 4), height = 30.0, twist = pi / 2)
mesh_slabs, = create_geometry(slabs)
plot_mesh(mesh_slabs; decorations = false)

# ## Disks and cylinders
#
# Disks in 2D ([`DiskSetup`](@ref)) with radii between `rmin` and `rmax`, possibly with
# several layers (here an axon and a myelin sheath), and an ECS: none ([`NoECS`](@ref)), a
# box ([`BoxECS`](@ref)) or the convex hull of the cells ([`ConvexHullECS`](@ref)), both
# with a margin.

disks = DiskSetup(;
    ncell = 6,
    rmin = 2.0,
    rmax = 5.0,
    nsidewall = 20,
    layersizes = [0.6, 1.0],
    ecs = BoxECS(; margin = 2.0),
)
mesh_disks, = create_geometry(disks; rng = Xoshiro(1))
plot_mesh(mesh_disks; colors = [fill(axon, 6); fill(myelin, 6); :lightgray])

# Cylinders ([`CylinderSetup`](@ref)) are disks extruded in 3D, and may be bent and
# twisted.

cylinders = CylinderSetup(;
    ncell = 4,
    rmin = 2.0,
    rmax = 4.0,
    nsidewall = 20,
    ecs = ConvexHullECS(; margin = 1.5),
    height = 30.0,
    bend = 0.02,
    twist = pi / 3,
)
mesh_cylinders, = create_geometry(cylinders; rng = Xoshiro(2))
plot_mesh(mesh_cylinders; colors = [fill(axon, 4); :gray], alpha = [ones(4); 0.25])

# ## Spheres
#
# Spheres ([`SphereSetup`](@ref)) may also have layers, for example a nucleus.

spheres = SphereSetup(;
    ncell = 4,
    rmin = 3.0,
    rmax = 5.0,
    layersizes = [0.5, 1.0],
    ecs = ConvexHullECS(; margin = 1.0),
)
mesh_spheres, = create_geometry(spheres; rng = Xoshiro(3))
colors = [fill(:crimson, 4); fill(:steelblue, 4); :gray]
plot_mesh(mesh_spheres; colors, alpha = [fill(1.0, 4); fill(0.4, 4); 0.15])

# ## Samples of dense tissue
#
# For dense tissue, the cells fill a periodic box ([`PeriodicBoxECS`](@ref)): the cells
# are packed randomly with gaps between them ([`periodic_packing`](@ref)), cut by the faces
# of the box, and their periodic images enter on the opposite faces. Here, 20 disks fill
# 70 % of the area.

packed_disks = DiskSetup(;
    ncell = 20,
    rmin = 1.0,
    rmax = 2.5,
    nsidewall = 20,
    dmin = 0.15,
    ecs = PeriodicBoxECS(; fraction = 0.7),
)
mesh_packed_disks, = create_geometry(packed_disks; rng = Xoshiro(4))
plot_mesh(mesh_packed_disks; colors = [fill(axon, 20); :lightgray])

# The same for 15 myelinated axons in 3D:

packed = CylinderSetup(;
    ncell = 15,
    rmin = 1.0,
    rmax = 2.5,
    nsidewall = 16,
    dmin = 0.15,
    layersizes = [0.7, 1.0],
    ecs = PeriodicBoxECS(; fraction = 0.65),
    height = 10.0,
)
mesh_packed, = create_geometry(packed; rng = Xoshiro(4))
plot_mesh(mesh_packed; colors = [fill(axon, 15); fill(myelin, 15); :lightgray])

# Polyhedral cells ([`VoronoiSetup`](@ref)) are the shrunk Voronoi cells of random seeds,
# and fill any fraction of the volume, here 85 %.

voronoi = VoronoiSetup(; ncell = 12, boxsize = 12.0, fraction = 0.85)
mesh_voronoi, = create_geometry(voronoi; rng = Xoshiro(5))
plot_mesh(mesh_voronoi; colors = [fill(:seagreen, 12); :lightgray], decorations = false)

# The meshes of these samples are periodic, for periodic boundary conditions (see
# [Dense tissue with periodic boundary conditions](@ref periodic)).
#
# Realistic neurons ([`NeuronSetup`](@ref)) are read from volume meshes, see
# [Neurons](@ref neuron).

# ## Overview
#
# All the recipes together:

fig = Figure(; size = (1200, 800))
panels = [
    (mesh_plates, "Plates", nothing),
    (
        mesh_disks,
        "Disks with layers and a box ECS",
        [fill(axon, 6); fill(myelin, 6); :lightgray],
    ),
    (mesh_slabs, "Twisted slabs", nothing),
    (mesh_cylinders, "Bent and twisted cylinders", [fill(axon, 4); (:gray, 0.25)]),
    (
        mesh_spheres,
        "Spheres with nuclei",
        [fill(:crimson, 4); fill((:steelblue, 0.4), 4); (:gray, 0.15)],
    ),
    (mesh_packed_disks, "Packed disks", [fill(axon, 20); :lightgray]),
    (
        mesh_packed,
        "Packed myelinated axons",
        [fill(axon, 15); fill(myelin, 15); :lightgray],
    ),
    (mesh_voronoi, "Polyhedral cells", [fill(:seagreen, 12); :lightgray]),
]
for (k, (mesh, title, colors)) in enumerate(panels)
    position = fig[div(k-1, 4)+1, mod1(k, 4)]
    dim = size(mesh.points[1], 1)
    ax =
        dim == 2 ? Axis(position; title, aspect = DataAspect()) :
        Axis3(position; title, aspect = :data)
    hidedecorations!(ax)
    hidespines!(ax)
    isnothing(colors) ? plot_mesh!(ax, mesh) : plot_mesh!(ax, mesh; colors)
end
fig
