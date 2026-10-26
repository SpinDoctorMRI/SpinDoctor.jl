# Apparent diffusion coefficients (ADC)

"""
    compute_adc_sta(model, gradient)

Compute the ADC of each compartment in the short diffusion time approximation (STA), for
the direction of the scalar gradient `gradient`.
"""
function compute_adc_sta(model::Model{T}, gradient::ScalarGradient) where {T}
    (; mesh, D) = model
    d = gradient.dir
    t = diffusion_time_sta(gradient.profile)

    adcs = zeros(T, length(mesh.points))
    for icmpt in eachindex(mesh.points)
        points = mesh.points[icmpt]
        facets = reduce(hcat, mesh.facets[icmpt, :])
        volume = sum(first(element_geometry(points, mesh.elements[icmpt])))
        areas, _, normals = facet_geometry(points, facets)
        D0 = d' * D[icmpt] * d

        # Surface area projected onto the plane orthogonal to the gradient direction
        area = sum((d' * normals)' .^ 2 .* areas)

        adcs[icmpt] = D0 * (1 - 4 / 3 * sqrt(D0 / pi) * area / volume * sqrt(t))
    end

    all(>=(0), adcs) ||
        @warn "Obtained negative STA ADC, STA approximation does not hold" filter(
            <(0),
            adcs,
        )

    adcs
end

"""
    solve_hadc(
        model,
        matrices,
        gradient;
        timestep = echotime(gradient) / 1000,
        arraytype = Array,
    )

Compute the ADC of each compartment in the direction of the scalar gradient `gradient`,
using the homogenized ADC model (HADC). See [`solve_hadc_tensors`](@ref).
"""
function solve_hadc(
    model,
    matrices,
    gradient::ScalarGradient;
    timestep = echotime(gradient) / 1000,
    arraytype = Array,
)
    tensors = solve_hadc_tensors(model, matrices, gradient.profile; timestep, arraytype)
    [gradient.dir' * D * gradient.dir for D in tensors]
end

"""
    solve_hadc_tensors(
        model,
        matrices,
        profile;
        timestep = echotime(profile) / 1000,
        arraytype = Array,
    )

Compute the apparent diffusion tensor of each compartment for the time profile `profile`,
using the homogenized ADC model (HADC) [Schiavi2016](@cite). The ADC in the direction `u` is
`u' * D * u`.

The HADC model assumes impermeable compartment boundaries: permeabilities, initial spin
densities and T2 relaxation are not taken into account. The ADC of a sample may be
estimated by weighting the compartment ADCs with [`compartment_weights`](@ref).

The auxiliary diffusion problem is solved with the Crank-Nicolson scheme, with steps of at
most `timestep` in each characteristic interval of the profile. By linearity, the ADC in any
direction follows from one solve per axis. With `arraytype = CuArray` (after `using CUDA,
CUDSS`), it is solved on an NVIDIA GPU.

With periodic matrices (`assemble_matrices(model; periodic = true)`), the compartments
are periodic, and the faces of the box are not boundaries.
"""
function solve_hadc_tensors(
    model,
    matrices,
    profile::AbstractTimeProfile;
    timestep = echotime(profile) / 1000,
    arraytype = Array,
)
    (; D) = model
    (; M_cmpts, S_cmpts, G, volumes) = matrices
    P_cmpts = get(matrices, :P_cmpts, nothing)
    if !isnothing(P_cmpts)
        # Periodic compartments: the matrices of the periodic points, where the surface
        # integrals on opposite faces of the box cancel
        M_cmpts = [sparse(P') * Mi * P for (Mi, P) in zip(M_cmpts, P_cmpts)]
        S_cmpts = [sparse(P') * Si * P for (Si, P) in zip(S_cmpts, P_cmpts)]
        G = [P' * Gi for (Gi, P) in zip(G, P_cmpts)]
    end

    # The auxiliary function solves M w' + S w = F(t) B u for the direction u, and the ADC
    # correction is u' A u / (volume * int_F2) with A = B' int F(t) w(t) dt. The compartments
    # are independent, and solved together (with block diagonal matrices).
    B = [G[i] * D[i] for i in eachindex(D)]
    Wint = hadc_integral(
        blockdiag(M_cmpts...),
        blockdiag(S_cmpts...),
        reduce(vcat, B),
        profile,
        timestep,
        arraytype,
    )
    offsets = cumsum([0; size.(B, 1)])
    map(eachindex(D)) do i
        A = B[i]' * Wint[(offsets[i]+1):offsets[i+1], :]
        D[i] - (A + A') / 2 / (volumes[i] * int_F2(profile))
    end
end

# Compute int_0^TE F(t) W(t) dt, where W solves M W' + S W = F(t) B, W(0) = 0, with the
# Crank-Nicolson scheme and the trapezoidal rule, on the device of `arraytype`
function hadc_integral(M, S, B, profile, timestep, arraytype)
    B = to_device(arraytype, B)
    W = zero(B)
    Wint = zero(B)
    rhs = similar(B)
    factorizations = Dict{eltype(B),Any}()
    ivals = intervals(profile)
    for k = 1:(length(ivals)-1)
        nstep = max(1, ceil(Int, (ivals[k+1] - ivals[k]) / timestep))
        h = (ivals[k+1] - ivals[k]) / nstep
        F = get!(() -> device_cholesky(arraytype, M + h / 2 * S), factorizations, h)
        E = to_device(arraytype, M - h / 2 * S)
        for n = 1:nstep
            t0 = ivals[k] + (n - 1) * h
            t1 = n == nstep ? ivals[k+1] : t0 + h
            f0 = integral(profile, t0)
            f1 = integral(profile, t1)
            Wint .+= h / 2 * f0 .* W
            mul!(rhs, E, W)
            rhs .+= h / 2 * (f0 + f1) .* B
            W = F \ rhs
            Wint .+= h / 2 * f1 .* W
        end
    end
    to_cpu(Wint)
end
