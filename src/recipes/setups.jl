# Setup recipes: geometries of cells and extracellular space

"""
    AbstractECS

Extracellular space (ECS) around the cells of a setup recipe.
"""
abstract type AbstractECS end

"""
    NoECS()

No extracellular space: the cells are the only compartments.
"""
struct NoECS <: AbstractECS end

"""
    BoxECS(; margin = 1.0)

Extracellular space filling the axis-aligned bounding box of the cells, extended by
`margin` on each side. The cells do not touch the faces of the box, so the ECS contains a
free layer of thickness `margin` along the faces.
"""
@with_kw struct BoxECS{T} <: AbstractECS
    margin::T = 1.0
    @assert margin > 0
end

"""
    ConvexHullECS(; margin = 1.0, threshold = pi / 20, maxiter = 10)

Extracellular space filling the convex hull of the cells, moved outwards by `margin`.
Facets of the hull with angles below `threshold` are refined for up to `maxiter` iterations
(3D only).
"""
@with_kw struct ConvexHullECS{T} <: AbstractECS
    margin::T = 1.0
    @assert margin > 0
    threshold::T = pi / 20
    @assert threshold >= 0
    maxiter::Int = 10
    @assert maxiter >= 0
end

"""
    PeriodicBoxECS(; fraction = 0.5)

Extracellular space of a sample of a periodic tissue. The box of the cell configuration
(`cells.box`, e.g. from [`periodic_packing`](@ref)) is filled with the cells, cut by its
faces, and their periodic images, and the ECS fills the rest. Unlike [`BoxECS`](@ref),
there is no free layer of ECS along the faces, and the geometry is periodic. The mesh is
periodic too (the opposite faces of the box have the same points), so that periodic
boundary conditions may be used with `assemble_matrices(model; periodic = true)`;
otherwise, the box faces are impermeable boundaries.

If no cell configuration is given to [`create_geometry`](@ref), a packing with the cell
volume fraction `fraction` is created (see [`create_cells`](@ref)). Available for disks
and cylinders. Fractions up to about 0.8 (counting the gaps as part of the cells) require
many cells, e.g. twenty or more. As in the other recipes, the disks are polygons inscribed in the circles
(with `nsidewall` sides for the mean radius), so that their area is slightly smaller than
that of the circles (by 1.1 % for 24 sides).
"""
@with_kw struct PeriodicBoxECS{T} <: AbstractECS
    fraction::T = 0.5
    @assert 0 < fraction < 1
end

"""
    AbstractSetup{T}

Setup recipe describing a geometry, from which [`create_geometry`](@ref) creates cells,
surfaces and a finite element mesh. All lengths are in µm.
"""
abstract type AbstractSetup{T} end

"""
    PlateSetup(; depth = 1.0, widths = [1.0, 1.0, 1.0, 1.0, 1.0], refinement = Inf)

Setup recipe for stacked plates in 2D, with the given `widths` along the x-axis and extent
`depth` along the y-axis. The triangles have an area of at most `refinement`.
"""
@with_kw struct PlateSetup{T} <: AbstractSetup{T}
    depth::T = 1.0
    @assert depth > 0
    widths::Vector{T} = [1.0, 1.0, 1.0, 1.0, 1.0]
    @assert all(>(0), widths)
    refinement::T = Inf
    @assert refinement > 0
end

"""
    DiskSetup(;
        ncell = 3,
        layersizes = [1.0],
        rmin = 2.0,
        rmax = 6.0,
        dmin = 0.2,
        dmax = 0.6,
        nsidewall = 12,
        ecs = NoECS(),
        refinement = Inf,
    )

Setup recipe for `ncell` disks in 2D, possibly immersed in an extracellular space `ecs`.

- The radii are between `rmin` and `rmax`, and the gaps between neighboring cells between
  `dmin` and `dmax` times the mean radius (see [`create_cells`](@ref)).
- Each cell has layers with the relative radii `layersizes` (increasing, the last one
  being one), e.g. `[0.6, 1.0]` for an axon with a myelin sheath.
- A cell of mean radius has `nsidewall` points on its boundary.
- The triangles have an area of at most `refinement`.
"""
@with_kw struct DiskSetup{T,E} <: AbstractSetup{T}
    ncell::Int = 3
    @assert ncell > 0
    layersizes::Vector{T} = [1.0]
    @assert all(s -> 0 < s <= 1, layersizes)
    @assert issorted(layersizes)
    rmin::T = 2.0
    rmax::T = 6.0
    @assert 0 < rmin <= rmax
    dmin::T = 0.2
    dmax::T = 0.6
    @assert 0 < dmin <= dmax
    nsidewall::Int = 12
    @assert nsidewall > 0
    ecs::E = NoECS()
    refinement::T = Inf
    @assert refinement > 0
end

"""
    SphereSetup(;
        ncell = 3,
        layersizes = [1.0],
        rmin = 2.0,
        rmax = 6.0,
        dmin = 0.2,
        dmax = 0.6,
        nsidewall = 200,
        ecs = NoECS(),
        refinement = Inf,
    )

Setup recipe for `ncell` spheres, possibly immersed in an extracellular space `ecs`. The
parameters are those of [`DiskSetup`](@ref), with `nsidewall` points on the surface of a
cell of mean radius, and a maximum tetrahedron volume `refinement`.
"""
@with_kw struct SphereSetup{T,E} <: AbstractSetup{T}
    ncell::Int = 3
    @assert ncell > 0
    layersizes::Vector{T} = [1.0]
    @assert all(s -> 0 < s <= 1, layersizes)
    @assert issorted(layersizes)
    rmin::T = 2.0
    rmax::T = 6.0
    @assert 0 < rmin <= rmax
    dmin::T = 0.2
    dmax::T = 0.6
    @assert 0 < dmin <= dmax
    nsidewall::Int = 200
    @assert nsidewall > 0
    ecs::E = NoECS()
    refinement::T = Inf
    @assert refinement > 0
end

"""
    NeuronSetup(; meshfile, ecs = NoECS(), refinement = Inf)

Setup recipe for a neuron, possibly wrapped in a ECS. The neuron is given by the volume
mesh `meshfile` in the Gmsh format, e.g. from the RealNeuronMeshes repository.
"""
@with_kw struct NeuronSetup{T,E} <: AbstractSetup{T}
    meshfile::String = ""
    ecs::E = NoECS()
    refinement::T = Inf
    @assert refinement > 0
end

"""
    ExtrusionSetup(; groundsetup, height = 1.0, bend = 0.0, twist = 0.0, refinement = Inf)

Setup recipe for the extrusion of the 2D setup `groundsetup` along the z-axis over `height`,
possibly bent and twisted (see [`deform_domain!`](@ref)). The tetrahedra have a volume of at
most `refinement`.
"""
@with_kw struct ExtrusionSetup{T,S} <: AbstractSetup{T}
    groundsetup::S
    refinement::T = Inf
    @assert refinement > 0
    height::T = 1.0
    @assert height > 0
    bend::T = 0.0
    twist::T = 0.0
end

"""
    CylinderSetup(;
        ncell = 1,
        layersizes = [1.0],
        rmin = 2.0,
        rmax = 6.0,
        dmin = 0.2,
        dmax = 0.6,
        nsidewall = 12,
        ecs = NoECS(),
        height = 1,
        bend = 0.0,
        twist = 0.0,
        refinement = Inf,
    )

Setup recipe for parallel cylinders: the extrusion of a [`DiskSetup`](@ref).
"""
CylinderSetup{T,E} = ExtrusionSetup{T,DiskSetup{T,E}}
CylinderSetup(;
    ncell = 1,
    layersizes = [1.0],
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.6,
    nsidewall = 12,
    ecs = NoECS(),
    height = 1.0,
    bend = 0.0,
    twist = 0.0,
    refinement = Inf,
) = ExtrusionSetup(;
    groundsetup = DiskSetup(; ncell, layersizes, rmin, rmax, dmin, dmax, nsidewall, ecs),
    height,
    bend,
    twist,
    refinement,
)

"""
    SlabSetup(;
        widths = [1.0, 1.0, 1.0, 1.0, 1.0],
        depth = 5.0,
        height = 1.0,
        bend = 0.0,
        twist = 0.0,
        refinement = Inf,
    )

Setup recipe for stacked slabs: the extrusion of a [`PlateSetup`](@ref).
"""
SlabSetup{T} = ExtrusionSetup{T,PlateSetup{T}}
SlabSetup(;
    widths = [1.0, 1.0, 1.0, 1.0, 1.0],
    depth = 5.0,
    height = 1.0,
    bend = 0.0,
    twist = 0.0,
    refinement = Inf,
) = ExtrusionSetup(;
    groundsetup = PlateSetup(; widths, depth),
    height,
    bend,
    twist,
    refinement,
)

"""
    radial_dimension(setup)

Dimension of the radial direction of plates (1), disks (2) and spheres (3).
"""
function radial_dimension end

radial_dimension(::PlateSetup) = 1
radial_dimension(::DiskSetup) = 2
radial_dimension(::SphereSetup) = 3

# Boundary numbering for n regions: the interfaces between regions i < j come first, row by
# row (12, 13, ..., 1n, 23, ...), followed by the outer boundaries of the regions 1, ..., n
ninterface(n) = div(n * (n - 1), 2)
nboundary(n) = ninterface(n) + n
interface_index(i, j, n) = div((i - 1) * (2n - i), 2) + j - i
outer_boundary_index(i, n) = ninterface(n) + i
