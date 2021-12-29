"""
    coefficients(setup; D, T₂, ρ, κ, γ)

Prepare PDE compartments.
"""
function coefficients(setup::AbstractSetup{T}; D, T₂, ρ, κ, γ) where T
    (; ncell, include_in, in_ratio, ecs_shape, ecs_ratio) = setup

    include_ecs = ecs_shape != "no_ecs"

    # Check for correct radius ratios and that neurons do not have in-compartments
    @assert !include_in || 0 < in_ratio && in_ratio < 1 && !isa(setup, NeuronSetup)
    @assert !include_ecs || 0 < ecs_ratio

    # Determine number of compartments and boundaries
    ncompartment = (1 + include_in) * ncell + include_ecs
    if isa(setup, SphereSetup)
        # For a sphere, there is one interface
        nboundary = (include_in + 1) * ncell + include_ecs
    elseif isa(setup, CylinderSetup)
        # An axon has a side interface, and a top-bottom boundary
        nboundary = (2 * include_in + 1 + include_ecs) * ncell + include_ecs
    elseif isa(setup, NeuronSetup)
        # For a neuron, there is one interface
        nboundary = 1 + include_ecs
    end

    boundaries = String[]
    boundary_markers = fill(false, ncompartment, nboundary)
    if include_in
        # Add in-compartments and in-out-interfaces
        compartments = repeat(["in"], ncell)
        boundaries = repeat(["in,out"], ncell)
        boundary_markers[CartesianIndex.(1:2*ncell, repeat(1:ncell, 2))] .= true
        ncompartment_old = ncell
        nboundary_old = ncell
    else
        compartments = []
        boundaries = []
        ncompartment_old = 0
        nboundary_old = 0
    end

    # Add out-compartments
    append!(compartments, repeat(["out"], ncell))

    if include_ecs
        # Add ecs-compartment and out-ecs interfaces
        push!(compartments, "ecs")
        append!(boundaries, repeat(["out,ecs"], ncell))
        boundary_markers[CartesianIndex.(
            [ncompartment_old+1:ncompartment_old+ncell; fill(ncompartment, ncell)],
            repeat(nboundary_old+1:nboundary_old+ncell, 2),
        )] .= true
        nboundary_old = nboundary_old + ncell
    end

    if isa(setup, CylinderSetup)
        if include_in
            # Add in boundary
            append!(boundaries, repeat(["in"], ncell))
            ncompartment_old = ncell
            boundary_markers[CartesianIndex.(
                1:ncell,
                nboundary_old+1:nboundary_old+ncell,
            )] .= true
            nboundary_old = nboundary_old + ncell
        else
            ncompartment_old = 0
        end

        # Add out boundary
        append!(boundaries, repeat(["out"], ncell))
        boundary_markers[CartesianIndex.(
            ncompartment_old+1:ncompartment_old+ncell,
            nboundary_old+1:nboundary_old+ncell,
        )] .= true

        if include_ecs
            # Add ecs boundary
            push!(boundaries, "ecs")
            boundary_markers[end, end] = true
        end
    elseif include_ecs
        # Add ecs boundary
        push!(boundaries, "ecs")
        boundary_markers[end, end] = true
    else
        # Add out boundary
        push!(boundaries, "out")
        boundary_markers[end, end] = true
    end

    # Initialize output arrays
    coeffs = (;
        D = [zeros(SMatrix{3,3,T}) for _ = 1:ncompartment],
        T₂ = zeros(T, ncompartment),
        κ = zeros(T, nboundary),
        ρ = zeros(Complex{T}, ncompartment),
        γ,
    ) 

    # Distribute material properties to compartments and boundaries
    coeffs.ρ[compartments.=="in"] .= ρ.in
    coeffs.ρ[compartments.=="out"] .= ρ.out
    coeffs.ρ[compartments.=="ecs"] .= ρ.ecs
    coeffs.D[compartments.=="in"] .= [D.in]
    coeffs.D[compartments.=="out"] .= [D.out]
    coeffs.D[compartments.=="ecs"] .= [D.ecs]
    coeffs.T₂[compartments.=="in"] .= T₂.in
    coeffs.T₂[compartments.=="out"] .= T₂.out
    coeffs.T₂[compartments.=="ecs"] .= T₂.ecs
    coeffs.κ[boundaries.=="in,out"] .= κ.in_out
    coeffs.κ[boundaries.=="out,ecs"] .= κ.out_ecs
    coeffs.κ[boundaries.=="in"] .= κ.in
    coeffs.κ[boundaries.=="out"] .= κ.out
    coeffs.κ[boundaries.=="ecs"] .= κ.ecs

    coeffs 
end
