"""
    mesh, surfaces, cells = create_geometry(
        setup;
        cells = nothing,
        rng = Random.default_rng(),
        cache = nothing,
        savedir = nothing,
        recreate = false,
    )

Create the cells, the surface triangulation and the finite element mesh of `setup`.

Random cell configurations are drawn from `rng`. For disks, cylinders and spheres, a cell
configuration `(; centers, radii)` (as returned by `create_cells`) may be given instead,
with `centers` of size `(dim, ncell)` and `radii` of length `ncell`, e.g. from a packing
algorithm.

By default, no files are written. The geometry may be stored in one of two ways:

- `cache`: a JLD2 file (requires `using JLD2`, see [`save_geometry`](@ref)). If the file
  exists and was created for the same setup (and cells), the geometry is loaded from it.
  Otherwise, it is created and saved.
- `savedir`: a directory with text files in the TetGen format, compatible with the MATLAB
  version of SpinDoctor (`cells`, `mesh.node`, `mesh.poly`, `mesh.1.node`, `mesh.1.face`,
  `mesh.1.ele`), and a description of the setup (`setup`). Existing files are reused if
  they were created for the same setup.

With `recreate = true`, the geometry is always created, and the files are overwritten.
"""
function create_geometry(
    setup;
    cells = nothing,
    rng = Random.default_rng(),
    cache = nothing,
    savedir = nothing,
    recreate = false,
)
    isnothing(cache) || isnothing(savedir) || error("Give either `cache` or `savedir`")
    if !isnothing(cache)
        if !recreate && isfile(cache)
            g = load_geometry(cache)
            if g.setup == repr(setup) && (isnothing(cells) || g.cells == cells)
                return g.mesh, g.surfaces, g.cells
            end
            @info "The geometry in $cache was created for a different setup. Recreating it."
        end
        mesh, surfaces, cells = create_geometry(setup; cells, rng)
        save_geometry(cache, setup, mesh, surfaces, cells)
        return mesh, surfaces, cells
    end
    isnothing(savedir) || return create_geometry_text(setup, savedir; cells, rng, recreate)

    cells = setup_cells(setup, cells, rng)
    if setup isa NeuronSetup
        isempty(setup.meshfile) && error("Give the volume mesh of the neuron as `meshfile`")
        surfaces = create_surfaces(setup, setup.meshfile)
    else
        surfaces = create_surfaces(setup, cells)
    end
    mesh_all = create_mesh(surfaces, setup.refinement)
    finalize_geometry(setup, mesh_all, surfaces, cells)
end

# Cell configuration: given or random
function setup_cells(setup, cells, rng)
    isnothing(cells) && return create_cells(setup; rng)
    ncell = cells_setup(setup).ncell
    size(cells.centers, 2) == length(cells.radii) == ncell ||
        error("The setup has $ncell cells, but $(length(cells.radii)) were given")
    cells
end

# Create the geometry, or reuse it, from text files in the TetGen format
function create_geometry_text(setup, savedir; cells, rng, recreate)
    mkpath(savedir)
    setupfile = joinpath(savedir, "setup")
    cellfile = joinpath(savedir, "cells")
    filename = joinpath(savedir, "mesh")

    reuse = !recreate && isnothing(cells) && isfile(filename * ".1.node")
    if reuse && isfile(setupfile) && read(setupfile, String) != repr(setup)
        @info "The files in $savedir were created for a different setup. Recreating them."
        reuse = false
    end

    if reuse
        cells = isfile(cellfile) ? read_cells(cellfile) : nothing
        surfaces = read_surfaces(filename)
        mesh_all = read_mesh(filename * ".1")
    else
        cells = setup_cells(setup, cells, rng)
        isnothing(cells) || save_cells(cellfile, cells)
        if setup isa NeuronSetup
            meshfile = isempty(setup.meshfile) ? filename * ".msh" : setup.meshfile
            surfaces = create_surfaces(setup, meshfile)
        else
            surfaces = create_surfaces(setup, cells)
        end
        save_surfaces(filename, surfaces)
        mesh_all = create_mesh(surfaces, setup.refinement)
        save_mesh(filename * ".1", mesh_all)
        write(setupfile, repr(setup))
    end

    finalize_geometry(setup, mesh_all, surfaces, cells)
end

# Deform the global mesh, and split it into compartments
function finalize_geometry(setup, mesh_all, surfaces, cells)
    if setup isa ExtrusionSetup && !(iszero(setup.bend) && iszero(setup.twist))
        @debug "Deforming domain with bend $(setup.bend) and twist $(setup.twist)"
        deform_domain!(mesh_all.points, setup.bend, setup.twist)
    end

    mesh = split_mesh(mesh_all)

    # Check that boundaries are correctly assigned
    ncompartment = length(mesh.points)
    for i = 1:ncompartment, j in [1:(i-1); (i+1):ncompartment]
        isempty(mesh.facets[i, end-ncompartment+j]) || @warn("""
             Outer boundary of compartment $i is registered as
             an interface to compartment $j, this is not intended.
             Consider rerunning `create_geometry`.
        """)
    end

    mesh, surfaces, cells
end
