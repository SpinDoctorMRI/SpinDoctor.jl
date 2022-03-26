"""
    create_geometry(setup; recreate = true)

Create cells, surfaces and finite element mesh. If `recreate = false`, previous geometry
will be reused.

This function does the following:
- Check geometry setup consistency
- Create or load cell configuration
- Create or load surface triangulation
- Call TetGen
- Deform domain
- Split mesh into compartments
For custom geometries with more than one compartment, call `split_mesh`
directly instead. This requires facet and element labels.
"""
function create_geometry(setup::Union{CylinderSetup,SphereSetup,NeuronSetup,EMSetup,FiberSetup,PlateSetup}; recreate = true)
    (; refinement, tetgen_option) = setup

    # File name for saving or loading geometry
    filename = joinpath(setup.meshdir, setup.name)

    # Make sure that folder exists
    dir = dirname(filename)
    isdir(dir) || mkpath(dir)

    # if hasfield(typeof(setup), :ecs_shape)
    #     setup.ecs_shape ∈ [:no_ecs, :box, :convex_hull, :tight_wrap] ||
    #         error("Invalid ECS shape")
    # end
    if isa(setup, PlateSetup)
        no_ecs = true
    else
        no_ecs = setup.ecs_shape == :no_ecs
    end

    if setup isa Union{CylinderSetup,SphereSetup}
        # Check if cell description file is already available
        cellfilename = filename * "_cells"
        if isfile(cellfilename) && !recreate
            cells = read_cells(cellfilename)
        else
            cells = create_cells(setup)
            save_cells(cells, cellfilename)
        end
    else
        cells = nothing
    end

    # Make directory for storing finite elements mesh
    is_stl = endswith(filename, ".stl")
    if is_stl
        setup.ecs_shape == :no_ecs || error("ECS is only available for surface meshes")
        filename = filename[1:end-4]
    end
    save_meshdir_path = filename * "_dir"
    isdir(save_meshdir_path) || mkpath(save_meshdir_path)

    # Use an existing finite elements mesh or create a new finite elements mesh. The name of
    # the finite elements mesh is stored in the string `fname_tetgen_femesh`
    refinement_str = isinf(refinement) ? "" : "_refinement$refinement"
    refinement_str = isempty(tetgen_option) ? refinement_str : refinement_str*"_tet$tetgen_option"
    ecs_str = no_ecs ? "" : "_$(setup.ecs_shape)$(setup.ecs_ratio)"
    fname_tetgen = "$save_meshdir_path/$(split(filename, "/")[end])$(ecs_str)$(refinement_str)_mesh"

    # Read or create surface triangulation
    if is_stl
        surfaces = nothing
    elseif isfile(fname_tetgen * ".node") && isfile(fname_tetgen * ".poly") && !recreate
        surfaces = read_surfaces(fname_tetgen)
    else
        if isa(setup, NeuronSetup)
            surfaces = create_surfaces(setup, filename)
        elseif isa(setup, EMSetup)
            surfaces = create_surfaces(setup, filename)
        elseif isa(setup, FiberSetup)
            surfaces = create_surfaces(setup, filename)
        else
            surfaces = create_surfaces(setup, cells)
        end
        save_surfaces(fname_tetgen, surfaces)
    end

    # Add ".1" suffix to output file name, since this is what Tetgen does
    fname_tetgen_femesh = fname_tetgen * ".1"

    if isfile(fname_tetgen_femesh * ".node") && !recreate
        # Read global mesh from Tetgen output
        mesh_all = read_tetgen(fname_tetgen_femesh)
        @assert size(surfaces.regions,2) == length(unique(mesh_all.elementmarkers)) "Incorrect number of compartments. $(length(surfaces.regions)) ≠ $(length(unique(mesh_all.elementmarkers))) "
        @assert length(unique(surfaces.facetmarkers)) == length(unique(mesh_all.facetmarkers)) "Incorrect number of boundaries. $(length(unique(surfaces.facetmarkers))) ≠ $(length(unique(mesh_all.facetmarkers)))"
    elseif is_stl
        error("Not implemented")
    else
        mesh_all = call_tetgen(surfaces, refinement, tetgen_option)
        save_tetgen(mesh_all, fname_tetgen_femesh)
    end

    # Deform domain
    if isa(setup, CylinderSetup) && (setup.deform_angle.bend > 1e-16 || setup.deform_angle.twist > 1e-16)
        @debug "Deforming domain with bend $(setup.bend) and twist $(setup.twist)"
        deform_domain!(mesh_all.points, setup.deform_angle.bend, setup.deform_angle.twist)
    end

    # Split mesh into compartments
    mesh = split_mesh(mesh_all)

    mesh, surfaces, cells
end
