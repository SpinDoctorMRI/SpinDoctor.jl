function initialize!(writer::VTKWriter, simulation, gradient, ξ, t)
    (; dir, filename) = writer 
    ispath(dir) || mkdir(dir)
    writer.pvd = paraview_collection(joinpath(dir, filename))
    update!(writer, simulation, gradient, ξ, t)
end

function initialize!(p::Plotter{T}, simulation, gradient, ξ, t) where {T}
    femesh = simulation.model.mesh
    M = simulation.matrices.M
    ncompartment, nboundary = size(femesh.facets)
    ξ_cmpts = split_field(femesh, ξ)
    ξmax = maximum(abs, ξ)
    p.S₀ = sum(abs, M * ξ)
    TE = echotime(gradient.profile)

    p.t[] = [t]
    p.f[] = [gradient.profile(t)]
    p.ξ[] = copy(ξ)
    p.attenuation = Node(T[1])
    p.magnitude = [Node(T[]) for _ = 1:ncompartment, _ = 1:nboundary]
    p.phase = [Node(T[]) for _ = 1:ncompartment, _ = 1:nboundary]

    fig = Figure()

    ax = Axis(fig[1, 1]; xlabel = "t [μs]", title = "Time profile")
    xlims!(ax, 0, TE)
    ylims!(ax, -1.1, 1.1)
    lines!(ax, p.t, p.f)

    # ax = Axis(fig[2, 1]; xlabel = "t [μs]", yscale = log10, title = "Signal")
    ax = Axis(fig[2, 1]; xlabel = "t [μs]", title = "Signal attenuation")
    xlims!(ax, 0, TE)
    # ylims!(ax, 1e-8, 1)
    ylims!(ax, 0, 1.1)
    lines!(ax, p.t, p.attenuation)

    gm = fig[1, 2] = GridLayout()
    ax = Axis3(gm[1, 1], title = "Magnetization (magnitude)")
    ax.aspect = :data 
    m = nothing
    first = true
    colorrange = (0, ξmax)
    for icmpt = 1:ncompartment, iboundary = 1:nboundary
        facets = femesh.facets[icmpt, iboundary]
        points = femesh.points[icmpt]
        color = p.magnitude[icmpt, iboundary]
        color[] = abs.(ξ_cmpts[icmpt])
        if first
            m = mesh!(ax, points', facets'; color, shading = false, colorrange)
            first = false
        else
            mesh!(ax, points', facets'; color, shading = false, colorrange)
        end
    end
    Colorbar(gm[1, 2]; limits = colorrange)

    gm = fig[2, 2] = GridLayout()
    ax = Axis3(gm[1, 1], title = "Magnetization (phase-shift)")
    ax.aspect = :data 
    m = nothing
    first = true
    # colormap = :cyclic_wrwbw_40_90_c42_n256
    colormap = :cyclic_mrybm_35_75_c68_n256		
    colorrange = (-π, π)
    for icmpt = 1:ncompartment, iboundary = 1:nboundary
        facets = femesh.facets[icmpt, iboundary]
        points = femesh.points[icmpt]
        color = p.phase[icmpt, iboundary]
        color[] = angle.(ξ_cmpts[icmpt])
        if first
            m = mesh!(ax, points', facets'; color, shading = false, colorrange, colormap)
            first = false
        else
            mesh!(ax, points', facets'; color, shading = false, colorrange, colormap)
        end
    end
    Colorbar(gm[1, 2]; limits = colorrange, colormap)

    display(fig)

    p
end
