"""
    analytical_coefficients(setup, coeffs)

Get coefficients for the analytical module.
"""
analytical_coefficients(s::ExtrusionSetup, coeffs) =
    analytical_coefficients(s.groundsetup, coeffs)

function analytical_coefficients(setup::PlateSetup, coeffs)
    widths = setup.widths
    dim = radial_dimension(setup)
    n = length(widths)

    ρ = coeffs.ρ
    r = cumsum(widths) * 1e-6
    D = tr.(coeffs.D) ./ size(coeffs.D[1], 1) .* 1e-6
    T₂ = coeffs.T₂ * 1e-6
    W = [[coeffs.κ[(i-1)*(2n-i)÷2+1] for i = 1:n-1]; coeffs.κ[end]]
    γ = coeffs.γ

    (; ρ, r, D, W, T₂, γ, dim)
end

function analytical_coefficients(setup::Union{DiskSetup,SphereSetup}, coeffs)
    layersizes = setup.layersizes
    dim = radial_dimension(setup)
    n = length(layersizes)

    r_mean = (setup.rmin + setup.rmax) / 2
    r = layersizes * rmean

    ρ = coeffs.ρ
    r = r * 1e-6
    D = tr.(coeffs.D) ./ size(coeffs.D[1], 1) .* 1e-6
    T₂ = coeffs.T₂ * 1e-6
    W = [[coeffs.κ[(i-1)*(2n-i)÷2+1] for i = 1:n-1]; coeffs.κ[end]]
    γ = coeffs.γ

    (; ρ, r, D, W, T₂, γ, dim)
end
