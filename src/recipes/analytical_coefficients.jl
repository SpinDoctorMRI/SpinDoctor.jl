function analytical_coefficients(setup, coeffs)
    rmean = (setup.rmin + setup.rmax) / 2
    include_in = setup.include_in
    include_ecs = setup.ecs_shape != "no_ecs"
    dim = radial_dimension(setup)
    n = unique([1:include_in+include_ecs; length(coeffs.κ)])

    # Get OUT parameters
    r_out = (setup.rmin + setup.rmax) / 2

    # Get IN parameters
    if setup.include_in
        r_in = setup.in_ratio * r_out
    else
        r_in = []
    end

    # Get ECS parameters
    if include_ecs 
        # Include ECS
        r_ecs = r_out + setup.ecs_ratio * rmean
    else
        # Empty ECS
        r_ecs = []
    end

    ρ = coeffs.ρ
    r = vcat(r_in, r_out, r_ecs) * 1e-6
    D = tr.(coeffs.D) ./ 3 .* 1e-6
    T₂ = coeffs.T₂ * 1e-6
    W = coeffs.κ[n]
    γ = coeffs.γ

    (; ρ, r, D, W, T₂, γ, dim)
end
