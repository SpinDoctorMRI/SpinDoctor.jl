"""
    create_geometry(setup; meshdir = nothing, recreate = true)

Create cells, surfaces and finite element mesh. If `meshdir` is a path,
geometry files are saved. If additionally `recreate = false`, previous
geometry files will be reused instead of generating new ones.
"""
function create_geometry(setup; meshdir = nothing, recreate = true)
    (; refinement) = setup

    # File name for saving or loading geometry
    if isnothing(meshdir)
        do_save = false
        meshdir = ""
    else
        do_save = true

        # Make sure that folder exists
        isdir(meshdir) || mkpath(meshdir)
    end

    if hasfield(typeof(setup), :ecs_shape)
        setup.ecs_shape ∈ [:no_ecs, :box, :convex_hull, :tight_wrap] ||
            error("Invalid ECS shape")
    end

    if setup isa Union{CylinderSetup,SphereSetup}
        # Check if cell description file is already available
        cellfilename = joinpath(meshdir, "cells")
        if do_save && isfile(cellfilename) && !recreate
            cells = read_cells(cellfilename)
        else
            cells = create_cells(setup)
            do_save && save_cells(cells, cellfilename)
        end
    else
        cells = nothing
    end

    stl_file = joinpath(meshdir, "mesh.stl")
    is_stl = do_save && isfile(stl_file)
    if is_stl
        setup.ecs_shape == :no_ecs || error("ECS is only available for surface meshes")
    end

    # Use an existing finite elements mesh or create a new finite elements mesh. The name of
    # the finite elements mesh is stored in the string `fname_tetgen_femesh`
    fname_tetgen = "$meshdir/mesh"

    # Read or create surface triangulation
    if is_stl
        surfaces = nothing
    elseif do_save &&
           isfile(fname_tetgen * ".node") &&
           isfile(fname_tetgen * ".poly") &&
           !recreate
        surfaces = read_surfaces(fname_tetgen)
    else
        if isa(setup, NeuronSetup)
            surfaces = create_surfaces(setup, fname_tetgen)
        else
            surfaces = create_surfaces(setup, cells)
        end
        do_save && save_surfaces(fname_tetgen, surfaces)
    end

    # Add ".1" suffix to output file name, since this is what Tetgen does
    fname_tetgen_femesh = fname_tetgen * ".1"

    if do_save && isfile(fname_tetgen_femesh * ".node") && !recreate
        # Read global mesh from Tetgen output
        mesh_all = read_mesh(fname_tetgen_femesh)
    elseif do_save && is_stl
        error("Not implemented")
    else
        mesh_all = create_mesh(surfaces, refinement)
        do_save && save_mesh(mesh_all, fname_tetgen_femesh)
    end

    # Deform domain
    if isa(setup, ExtrusionSetup) && !(setup.bend ≈ 0 && setup.twist ≈ 0)
        @debug "Deforming domain with bend $(setup.bend) and twist $(setup.twist)"
        deform_domain!(mesh_all.points, setup.bend, setup.twist)
    end

    # Split mesh into compartments
    mesh = split_mesh(mesh_all)

    ncompartment = length(mesh.points)
    for i = 1:ncompartment, j in [1:i-1; i+1:ncompartment]
        isempty(mesh.facets[i, end-ncompartment+j]) || @warn("""
             Outer boundary of compartment $i is registered as
             an interface to compartment $j, this is not intended.
             Consider rerunning `create_geometry`.
        """)
    end

    mesh, surfaces, cells
end
