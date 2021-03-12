"""
    create_geometry(setup)

Create cells, surfaces and finite element mesh.
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
function create_geometry(setup::Setup)

    @unpack ncell, cell_shape, ecs_shape, deformation = setup.geometry
    @unpack compartments, boundaries = setup.pde

    refinement = haskey(setup.geometry, :refinement) ? setup.geometry[:refinement] : nothing

    # File name for saving or loading geometry
    filename = "meshfiles" * "/" * setup.name

    # Make sure that folder exists
    dir = dirname(filename)
    isdir(dir) || mkdir(dir)

    # Check correct input format
    @assert cell_shape ∈ ["sphere", "cylinder", "neuron"]
    if cell_shape == "neuron"
        # We do not deform the finite elements mesh of neurons
        if haskey(setup.geometry, :deformation) && any(deformation .> 1e-16)
            @error "Deformation is not available for neurons. Set deformation to (0, 0)."
        end
        if ncell != 1
            error("Neuron cell type is only available for ncell=1.")
        end
    end
    !(ecs_shape == "no_ecs" && ncell > 1) ||
        @error "Geometry must include ECS if more than one cell"
    @assert ecs_shape ∈ ["no_ecs", "box", "convex_hull", "tight_wrap"]

    # Check if cell description file is already available
    cellfilename = filename * "_cells"
    if isfile(cellfilename)
        cells = read_cells(cellfilename)
    elseif cell_shape ∈ ["sphere", "cylinder"]
        # Create cells
        cells = create_cells(setup)
        
        # Save cell configuration
        save_cells(cells, cellfilename)
    else
        # No cells
        cells = nothing
    end


    # Make directory for storing finite elements mesh
    is_stl = endswith(filename, ".stl")
    if is_stl
        # ECS is currently only available for surface meshes
        @assert ecs_shape == "no_ecs"
        filename = filename[1:end-4]
    end
    save_meshdir_path = filename * "_dir"
    if !isdir(save_meshdir_path)
        mkdir(save_meshdir_path)
    end
    
    
    # Use an existing finite elements mesh or create a new finite
    # elements mesh. The name of the finite elements mesh is stored in the string
    # fname_tetgen_femesh
    refinement_str = isnothing(refinement) ? "" : "_refinement$(setup.geometry[:refinement])"
    fname_tetgen = save_meshdir_path * "/" * split(filename, "/")[end] * refinement_str * "_mesh"
    
    # Read or create surface triangulation
    if isfile(fname_tetgen * ".node") && isfile(fname_tetgen * ".poly")
        surfaces = read_surfaces(fname_tetgen)
    else
        if cell_shape ==  "sphere"
            # Create surface geometry of spheres
            surfaces = create_surfaces_sphere(cells, setup)
        elseif cell_shape ==  "cylinder"
            # Create surface geometry of cylinders
            surfaces = create_surfaces_cylinder(cells, setup)
        elseif cell_shape ==  "neuron"
            if is_stl
                surfaces = nothing
            else
                surfaces = create_surfaces_neuron(filename, setup)
            end
        end
    
        if !is_stl
            save_surfaces(fname_tetgen, surfaces);
        end
    end

    # Add ".1" suffix to output file name, since this is what Tetgen does
    fname_tetgen_femesh = fname_tetgen * ".1";
    
    if isfile(fname_tetgen_femesh * ".node")
        # Read global mesh from Tetgen output
        mesh_all = read_tetgen(fname_tetgen_femesh)
    elseif is_stl
        error("Not implemented")
    else
        mesh_all = call_tetgen(surfaces, refinement)
        save_tetgen(mesh_all, fname_tetgen_femesh)
    end
    
    # Check that at correct number of compartments and boundaries has been found
    compartments_new = unique(mesh_all.elementmarkers)
    boundaries_new = unique(mesh_all.facetmarkers)

    solution = "use smaller refinement or change surface triangulation."
    length(compartments_new) == length(compartments) || @error "Incorrect number of compartments, " * solution
    length(boundaries_new)  == length(boundaries) || @error "Incorrect number of boundaries, " * solution
    
    # Deform domain
    if any(deformation .> 1e-16)
        @printf "Deforming domain with bend %g and twist %g\n" deformation...
        deform_domain!(mesh_all.points, deformation)
    end
    
    # Split mesh into compartments
    mesh = split_mesh(mesh_all, setup)

    mesh, surfaces, cells
end