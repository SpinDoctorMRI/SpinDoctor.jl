function initialize!(p::Printer, problem, gradient, ξ, t)
    @info "Solving for" typeof(problem) gradient
end

function initialize!(writer::VTKWriter, problem, gradient, ξ, t)
    (; dir, filename) = writer
    ispath(dir) || mkdir(dir)
    writer.pvd = paraview_collection(joinpath(dir, filename))
    update!(writer, problem, gradient, ξ, t)
end

function initialize!(p::Plotter{T,dim}, problem, gradient, ξ, t) where {T,dim}
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
    p.attenuation = Observable(T[1])
    if dim == 2
        p.magnitude = [Observable(T[]) for _ = 1:ncompartment, _ = 1:1]
        p.phase = [Observable(T[]) for _ = 1:ncompartment, _ = 1:1]
    elseif dim == 3
        p.magnitude = [Observable(T[]) for _ = 1:ncompartment, _ = 1:nboundary]
        p.phase = [Observable(T[]) for _ = 1:ncompartment, _ = 1:nboundary]
    end

    p.fig = Figure()

    if gradient isa ScalarGradient
        p.f[] = [gradient.profile(t)]
        ax = Axis(p.fig[1, 1]; xlabel = "t [μs]", title = "Time profile")
        xlims!(ax, 0, TE)
        ylims!(ax, -1.1, 1.1)
        lines!(ax, p.t, p.f)
    else
        grads = mapreduce(gradient, hcat, LinRange(0, TE, 200))
        pmin = min.(minimum(grads; dims = 2), zero(T))
        pmax = max.(maximum(grads; dims = 2), zero(T))
        inds = pmin .≈ pmax
        gmax = maximum(norm, grads)
        pmin[inds] .-= gmax / 2
        pmax[inds] .+= gmax / 2
        grad = Vec{dim,Float32}(gradient(t))
         p.g⃗ = Observable([grad])
         p.g⃗_hist = Observable([grad])
        if dim == 2
            ax = Axis(p.fig[1, 1]; title = "Gradient [T/m]")
            ax.aspect = :data
            xlims!(ax, pmin[1], pmax[1])
            ylims!(ax, pmin[2], pmax[2])
             lines!(ax, p.g⃗_hist)
             arrows!(ax, [Point2f(0, 0)], p.g⃗)
        elseif dim == 3
            ax = Axis3(p.fig[1, 1]; title = "Gradient [T/m]")
            ax.aspect = :data
            xlims!(ax, pmin[1], pmax[1])
            ylims!(ax, pmin[2], pmax[2])
            zlims!(ax, pmin[3], pmax[3])
             lines!(ax, p.g⃗_hist)
             arrows!(ax, [Point3f(0, 0, 0)], p.g⃗)
        end
    end

    ax = Axis(p.fig[2, 1]; xlabel = "t [μs]", title = "Signal attenuation")
    xlims!(ax, 0, TE)
    # ylims!(ax, 1e-8, 1)
    ylims!(ax, 0, 1.1)
    lines!(ax, p.t, p.attenuation)

    if dim == 2
        gm = p.fig[1, 2] = GridLayout()
        ax = Axis(gm[1, 1]; title = "Magnetization (magnitude)")
        colorrange = (0, ξmax)
        for icmpt = 1:ncompartment
            elements = femesh.elements[icmpt]
            points = femesh.points[icmpt]
            color = p.magnitude[icmpt]
            color[] = abs.(ξ_cmpts[icmpt])
            mesh!(ax, points', elements'; color, shading = false, colorrange)
        end
    elseif dim == 3
        gm = p.fig[1, 2] = GridLayout()
        ax = Axis3(gm[1, 1]; title = "Magnetization (magnitude)")
        ax.aspect = :data
        colorrange = (0, ξmax)
        for icmpt = 1:ncompartment, iboundary = 1:nboundary
            facets = femesh.facets[icmpt, iboundary]
            if !isempty(facets)
                points = femesh.points[icmpt]
                color = p.magnitude[icmpt, iboundary]
                color[] = abs.(ξ_cmpts[icmpt])
                mesh!(ax, points', facets'; color, shading = false, colorrange)
            end
        end
    end
    Colorbar(gm[1, 2]; limits = colorrange)

    gm = p.fig[2, 2] = GridLayout()
    # colormap = :cyclic_wrwbw_40_90_c42_n256
    colormap = :cyclic_mrybm_35_75_c68_n256
    # colorrange = (-0.99π, 0.99π)
    colorrange = (-π, π)
    if dim == 2
        gm = p.fig[2, 2] = GridLayout()
        ax = Axis(gm[1, 1]; title = "Magnetization (phase-shift)")
        for icmpt = 1:ncompartment
            elements = femesh.elements[icmpt]
            points = femesh.points[icmpt]
            color = p.phase[icmpt]
            color[] = angle.(ξ_cmpts[icmpt])
            mesh!(ax, points', elements'; color, shading = false, colorrange, colormap)
        end
    elseif dim == 3
        ax = Axis3(gm[1, 1]; title = "Magnetization (phase-shift)")
        ax.aspect = :data
        for icmpt = 1:ncompartment, iboundary = 1:nboundary
            facets = femesh.facets[icmpt, iboundary]
            if !isempty(facets)
                points = femesh.points[icmpt]
                color = p.phase[icmpt, iboundary]
                color[] = angle.(ξ_cmpts[icmpt])
                mesh!(ax, points', facets'; color, shading = false, colorrange, colormap)
            end
        end
    end
    Colorbar(gm[1, 2]; limits = colorrange, colormap)

    display(p.fig)

    p
end
