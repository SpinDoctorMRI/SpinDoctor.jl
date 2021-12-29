function initialize!(p::Printer, problem, gradient, ξ, t)
    @info "Solving for" typeof(problem) gradient
end

function initialize!(writer::VTKWriter, problem, gradient, ξ, t)
    (; dir, filename) = writer
    ispath(dir) || mkdir(dir)
    writer.pvd = paraview_collection(joinpath(dir, filename))
    update!(writer, problem, gradient, ξ, t)
end

function initialize!(p::Plotter{T}, problem, gradient, ξ, t) where {T}
    femesh = problem.model.mesh
    M = problem.matrices.M
    ncompartment, nboundary = size(femesh.facets)
    ξ_cmpts = split_field(femesh, ξ)
    ξmax = maximum(abs, ξ)
    p.S₀ = sum(abs, M * ξ)
    TE = echotime(gradient)

    p.n = 1
    p.t[] = [t]
    p.ξ[] = copy(ξ)
    p.attenuation = Node(T[1])
    p.magnitude = [Node(T[]) for _ = 1:ncompartment, _ = 1:nboundary]
    p.phase = [Node(T[]) for _ = 1:ncompartment, _ = 1:nboundary]

    p.fig = Figure()

    if gradient isa ScalarGradient
        p.f[] = [gradient.profile(t)]
        ax = Axis(p.fig[1, 1]; xlabel = "t [μs]", title = "Time profile")
        xlims!(ax, 0, TE)
        ylims!(ax, -1.1, 1.1)
        lines!(ax, p.t, p.f)
    else
        grads = mapreduce(gradient, hcat, LinRange(0, TE, 200))
        pmin = minimum(grads; dims = 2)
        pmax = maximum(grads; dims = 2)
        grad = Vec3f(gradient(t))
        p.g⃗ = Node([grad])
        p.g⃗_hist = Node([grad])
        ax = Axis3(p.fig[1, 1]; title = "Gradient [T/m]")
        ax.aspect = :data
        xlims!(ax, pmin[1], pmax[1])
        ylims!(ax, pmin[2], pmax[2])
        zlims!(ax, pmin[3], pmax[3])
        lines!(ax, p.g⃗_hist)
        arrows!(ax, [Point3f(0,0,0)], p.g⃗)
    end

    ax = Axis(p.fig[2, 1]; xlabel = "t [μs]", title = "Signal attenuation")
    xlims!(ax, 0, TE)
    # ylims!(ax, 1e-8, 1)
    ylims!(ax, 0, 1.1)
    lines!(ax, p.t, p.attenuation)

    gm = p.fig[1, 2] = GridLayout()
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

    gm = p.fig[2, 2] = GridLayout()
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

    display(p.fig)

    p
end
