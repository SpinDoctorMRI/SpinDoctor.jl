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
    p.n = 1
    p.t = Observable(t)
    p.ξ = Observable(ξ)

    femesh = problem.model.mesh
    M = problem.matrices.M
    ncompartment, nboundary = size(femesh.facets)
    dim = size(femesh.points[1], 1)

    ξmax = maximum(abs, ξ)
    S₀ = sum(abs, M * ξ)
    TE = echotime(gradient)

    npoint_cmpts = size.(femesh.points, 2)
    inds = [0; cumsum(npoint_cmpts[1:end])]
    ξ_cmpts = [@lift($(p.ξ)[1+inds[i]:inds[i+1]]) for i = 1:ncompartment]

    p.fig = Figure()

    if gradient isa ScalarGradient
        _f = Point2f[]
        f = @lift push!(_f, Point2f($(p.t), gradient.profile($(p.t))))
        ax = Axis(p.fig[1, 1]; xlabel = "t [μs]", title = "Time profile")
        xlims!(ax, 0, TE)
        ylims!(ax, -1.1, 1.1)
        lines!(ax, f)
    else
        grads = mapreduce(gradient, hcat, LinRange(0, TE, 200))
        pmin = min.(minimum(grads; dims = 2), zero(T))
        pmax = max.(maximum(grads; dims = 2), zero(T))
        i = pmin .≈ pmax
        gmax = maximum(norm, grads)
        pmin[i] .-= gmax / 2
        pmax[i] .+= gmax / 2
        gvec = @lift Vec{dim,Float32}(gradient($(p.t)))
        _gvecs = Vec{dim,Float32}[]
        gvecs = @lift push!(_gvecs, $gvec)
        if dim == 2
            ax = Axis(p.fig[1, 1]; title = "Gradient [T/m]")
            ax.aspect = :data
            xlims!(ax, pmin[1], pmax[1])
            ylims!(ax, pmin[2], pmax[2])
            lines!(ax, gvecs)
            arrows!(ax, [Point2f(0, 0)], @lift([$gvec]))
        elseif dim == 3
            ax = Axis3(p.fig[1, 1]; title = "Gradient [T/m]")
            ax.aspect = :data
            xlims!(ax, pmin[1], pmax[1])
            ylims!(ax, pmin[2], pmax[2])
            zlims!(ax, pmin[3], pmax[3])
            lines!(ax, gvecs)
            arrows!(ax, [Point3f(0, 0, 0)], @lift([$gvec]))
        end
    end

    # Note: `p.ξ[]` must already be updated before `p.t`
    _attenuation = Point2f[]
    attenuation = @lift push!(_attenuation, Point2f($(p.t), sum(abs, M * p.ξ[]) / S₀))

    ax = Axis(p.fig[2, 1]; xlabel = "t [μs]", title = "Signal attenuation")
    xlims!(ax, 0, TE)
    # ylims!(ax, 1e-8, 1)
    ylims!(ax, 0, 1.1)
    lines!(ax, attenuation)

    gm = p.fig[1, 2] = GridLayout()
    if dim == 2
        ax = Axis(gm[1, 1]; title = "Magnetization (magnitude)")
        colorrange = (0, ξmax)
        for icmpt = 1:ncompartment
            elements = femesh.elements[icmpt]
            points = femesh.points[icmpt]
            color = @lift abs.($(ξ_cmpts[icmpt]))
            mesh!(ax, points', elements'; color, shading = false, colorrange)
        end
    elseif dim == 3
        ax = Axis3(gm[1, 1]; title = "Magnetization (magnitude)")
        ax.aspect = :data
        colorrange = (0, ξmax)
        for icmpt = 1:ncompartment
            color = @lift abs.($(ξ_cmpts[icmpt]))
            points = femesh.points[icmpt]
            for iboundary = 1:nboundary
                facets = femesh.facets[icmpt, iboundary]
                if !isempty(facets)
                    mesh!(ax, points', facets';
                        color,
                        shading = false,
                        colorrange,
                    )
                end
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
            color = @lift angle.($(ξ_cmpts[icmpt]))
            elements = femesh.elements[icmpt]
            points = femesh.points[icmpt]
            mesh!(ax, points', elements'; color, shading = false, colorrange, colormap)
        end
    elseif dim == 3
        ax = Axis3(gm[1, 1]; title = "Magnetization (phase-shift)")
        ax.aspect = :data
        for icmpt = 1:ncompartment
            color = @lift angle.($(ξ_cmpts[icmpt]))
            points = femesh.points[icmpt]
            for iboundary = 1:nboundary
                facets = femesh.facets[icmpt, iboundary]
                if !isempty(facets)
                    mesh!(ax, points', facets'; color, shading = false, colorrange, colormap)
                end
            end
        end
    end
    Colorbar(gm[1, 2]; limits = colorrange, colormap)

    display(p.fig)

    p
end
