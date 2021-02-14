""" Prepare PDE compartments """
function prepare_pde(cellsetup::CellSetup, domainsetup::DomainSetup)
    @unpack shape, ncell, include_nucleus, nucleus_radiusratio, include_ecs, ecs_shape = cellsetup

    # Determine number of compartments and boundaries
    ncompartment = (1 + include_nucleus) * ncell + include_ecs
    if shape == "sphere"
        nboundary = ncell
    elseif shape == "cylinder"
        nboundary = 2ncell
    elseif shape == "neuron"
        nboundary = 1
    end
    nboundary = (1 + include_nucleus) * nboundary + include_ecs

    # Label compartments
    compartments = String[]
    include_nucleus && append!(compartments, repeat(["in"], ncell))
    append!(compartments, repeat(["out"], ncell))
    include_ecs && push!(compartments, "ecs")

    boundaries = String[]
    boundary_markers = fill(false, ncompartment, nboundary)
    if include_nucleus
        if shape == "cylinder"
            append!(boundaries, repeat(["in"], ncell))
            append!(boundaries, repeat(["out"], ncell))
            boundary_markers[CartesianIndex.(1:2ncell, 1:2ncell)] .= true
            nboundary_previous = 2ncell
        else
            nboundary_previous = 0
        end
        append!(boundaries, repeat(["in-out"], ncell))
        boundary_markers[CartesianIndex.(1:2ncell, [
            nboundary_previous+1:nboundary_previous+ncell
            nboundary_previous+1:nboundary_previous+ncell
        ])] .= true
        ncompartment_previous = ncell
        nboundary_previous = nboundary_previous + ncell
    elseif shape == "cylinder"
        append!(boundaries, repeat(["out"], ncell))
        boundary_markers[CartesianIndex.(1:ncell, 1:ncell)] .= true
        ncompartment_previous = 0
        nboundary_previous = ncell
    else
        ncompartment_previous = 0
        nboundary_previous = 0
    end
    if include_ecs
        append!(boundaries, repeat(["out-ecs"], ncell))
        push!(boundaries, "ecs")
        boundary_markers[CartesianIndex.(
            ncompartment_previous+1:ncompartment_previous+ncell,
            nboundary_previous+1:nboundary_previous+ncell
        )] .= true
        boundary_markers[ncompartment, nboundary_previous+1:nboundary] .= true
    else
        append!(boundaries, repeat(["out"], ncell))
        boundary_markers[CartesianIndex.(
            ncompartment_previous+1:ncompartment,
            nboundary_previous+1:nboundary
        )] .= true
    end

    # Initialize output arrays
    diffusivity = zeros(ncompartment)
    relaxation = zeros(ncompartment)
    permeability = zeros(nboundary)
    initial_density = zeros(ncompartment)

    # Distribute material properties to compartments and boundaries
    diffusivity[compartments .== "in"] .= domainsetup.diffusivity_in
    diffusivity[compartments .== "out"] .= domainsetup.diffusivity_out
    diffusivity[compartments .== "ecs"] .= domainsetup.diffusivity_ecs
    relaxation[compartments .== "in"] .= domainsetup.relaxation_in
    relaxation[compartments .== "out"] .= domainsetup.relaxation_out
    relaxation[compartments .== "ecs"] .= domainsetup.relaxation_ecs
    permeability[boundaries .== "in-out"] .= domainsetup.permeability_in_out
    permeability[boundaries .== "out-ecs"] .= domainsetup.permeability_out_ecs
    initial_density[compartments .== "in"] .= domainsetup.initial_density_in
    initial_density[compartments .== "out"] .= domainsetup.initial_density_out
    initial_density[compartments .== "ecs"] .= domainsetup.initial_density_ecs

    # Return named tuple
    (; boundary_markers, compartments, boundaries, ncompartment, nboundary, diffusivity, relaxation, permeability, initial_density)
end
