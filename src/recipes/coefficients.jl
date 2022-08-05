"""
    coefficients(setup; D, T₂, ρ, κ, γ)

Order coefficients compartment arrays.
"""
function coefficients end

coefficients(::PlateSetup{T}; D, T₂, ρ, κ, γ) where {T} = (;
    D = T.(D),
    T₂ = T.(T₂),
    κ = [T.(κ.interfaces); T.(κ.boundaries)],
    ρ = Complex{T}.(ρ),
    γ = T(γ),
)

function coefficients(setup::NeuronSetup{T}; D, T₂, ρ, κ, γ) where {T}
    (; ecs_shape, ecs_ratio) = setup

    include_ecs = ecs_shape != :no_ecs

    @assert !include_ecs || 0 < ecs_ratio

    # Determine number of compartments and boundaries
    ncompartment = 1 + include_ecs
    nboundary = 1 + include_ecs
    if include_ecs
        compartments = ["neuron", "ecs"]
        boundaries = ["neuron,ecs", "ecs"]
    else
        compartments = ["neuron"]
        boundaries = ["neuron"]
    end

    # Initialize output arrays
    coeffs = (;
        D = [zeros(T, 3, 3) for _ = 1:ncompartment],
        T₂ = zeros(T, ncompartment),
        κ = zeros(T, nboundary),
        ρ = zeros(Complex{T}, ncompartment),
        γ = T(γ),
    )

    # Distribute material properties to compartments and boundaries
    coeffs.ρ[compartments.=="neuron"] .= ρ.neuron
    coeffs.ρ[compartments.=="ecs"] .= ρ.ecs
    coeffs.D[compartments.=="neuron"] .= [D.neuron]
    coeffs.D[compartments.=="ecs"] .= [D.ecs]
    coeffs.T₂[compartments.=="neuron"] .= T₂.neuron
    coeffs.T₂[compartments.=="ecs"] .= T₂.ecs
    coeffs.κ[boundaries.=="neuron,ecs"] .= κ.neuron_ecs
    coeffs.κ[boundaries.=="neuron"] .= κ.neuron
    coeffs.κ[boundaries.=="ecs"] .= κ.ecs

    coeffs
end

function coefficients(
    setup::Union{CylinderSetup{T},SphereSetup{T}};
    D,
    T₂,
    ρ,
    κ,
    γ,
) where {T}
    (; ncell, include_in, in_ratio, ecs_shape, ecs_ratio) = setup

    include_ecs = ecs_shape != :no_ecs

    # Check for correct radius ratios and that neurons do not have in-compartments
    @assert !include_in || 0 < in_ratio && in_ratio < 1
    @assert !include_ecs || 0 < ecs_ratio

    # Determine number of compartments and boundaries
    ncompartment = (1 + include_in) * ncell + include_ecs
    if isa(setup, SphereSetup)
        # For a sphere, there is one interface
        nboundary = (include_in + 1) * ncell + include_ecs
        dim = 3
    elseif isa(setup, CylinderSetup)
        # An axon has a side interface, and a top-bottom boundary
        nboundary = (2 * include_in + 1 + include_ecs) * ncell + include_ecs
        dim = 3
    end

    boundaries = String[]
    if include_in
        # Add in-compartments and in-out-interfaces
        compartments = repeat(["in"], ncell)
        boundaries = repeat(["in,out"], ncell)
    else
        compartments = []
        boundaries = []
    end

    # Add out-compartments
    append!(compartments, repeat(["out"], ncell))

    if include_ecs
        # Add ecs-compartment and out-ecs interfaces
        push!(compartments, "ecs")
        append!(boundaries, repeat(["out,ecs"], ncell))
    end

    if isa(setup, CylinderSetup)
        if include_in
            # Add in boundary
            append!(boundaries, repeat(["in"], ncell))
        end

        # Add out boundary
        append!(boundaries, repeat(["out"], ncell))

        if include_ecs
            # Add ecs boundary
            push!(boundaries, "ecs")
        end
    elseif include_ecs
        # Add ecs boundary
        push!(boundaries, "ecs")
    else
        # Add out boundary
        push!(boundaries, "out")
    end

    # Initialize output arrays
    coeffs = (;
        D = [zeros(T, dim, dim) for _ = 1:ncompartment],
        T₂ = zeros(T, ncompartment),
        κ = zeros(T, nboundary),
        ρ = zeros(Complex{T}, ncompartment),
        γ = T(γ),
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
