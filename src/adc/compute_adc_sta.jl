"""
    compute_adc_sta(model)

Compute the ADC for each compartment in the short diffusion time regime.
"""
function compute_adc_sta(model::Model{T}, gradient::ScalarGradient) where {T}
    ncompartment = length(model.mesh.points)

    adcs = zeros(T, ncompartment)
    for icmpt = 1:ncompartment
        points = model.mesh.points[icmpt]
        elements = model.mesh.elements[icmpt]
        facets = reduce(hcat, model.mesh.facets[icmpt, :])

        volumes, = get_mesh_volumes(points, elements)
        v = sum(volumes)
        _, areas, _, normals = get_mesh_surface(points, facets)
        d = gradient.dir
        D₀ = d' * model.D[icmpt] * d

        # Project surface areas onto plane orthogonal to gradient direction
        # (result is 1×1 matrix, extract scalar component)
        SAu = ((d'normals) .^ 2*areas)[1]

        t = diffusion_time_sta(gradient.profile)
        adcs[icmpt] = D₀ * (1 - 4 / 3 * √(D₀ / π) * SAu / v * √t)
    end

    all(≥(0), adcs) ||
        @warn "Obtained negative STA ADC, STA approximation does not hold" filter(
            <(0),
            adcs,
        )

    adcs
end
