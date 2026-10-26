function initialize_plotter!(p::Plotter, problem, gradient, u, t)
    T = real(eltype(u))
    p.n = 1
    tobs = Observable(t)
    uobs = Observable(u)

    femesh = problem.model.mesh
    M = problem.matrices.M
    dim = size(femesh.points[1], 1)

    umax = maximum(abs, u)
    S0 = abs(sum(M * u))
    TE = echotime(gradient)

    fig = Figure(; size = (1000, 700))

    if gradient isa ScalarGradient
        _f = Point2f[]
        f = @lift push!(_f, Point2f($tobs, gradient.profile($tobs)))
        ax = Axis(fig[1, 1]; xlabel = "t [µs]", title = "Time profile")
        xlims!(ax, 0, TE)
        ylims!(ax, -1.1, 1.1)
        lines!(ax, f)
    else
        grads = mapreduce(gradient, hcat, LinRange(0, TE, 200))
        pmin = min.(minimum(grads; dims = 2), zero(T))
        pmax = max.(maximum(grads; dims = 2), zero(T))
        i = isapprox.(pmin, pmax)
        gmax = maximum(norm, eachcol(grads))
        pmin[i] .-= gmax / 2
        pmax[i] .+= gmax / 2
        gvec = @lift Vec{dim,Float32}(gradient($tobs))
        _gvecs = Point{dim,Float32}[]
        gvecs = @lift push!(_gvecs, Point{dim,Float32}($gvec))
        if dim == 2
            ax = Axis(fig[1, 1]; title = "Gradient [T/m]")
            ax.aspect = DataAspect()
            xlims!(ax, pmin[1], pmax[1])
            ylims!(ax, pmin[2], pmax[2])
            lines!(ax, gvecs)
            arrows2d!(ax, [Point2f(0, 0)], @lift([$gvec]))
        elseif dim == 3
            ax = Axis3(fig[1, 1]; title = "Gradient [T/m]")
            ax.aspect = :data
            xlims!(ax, pmin[1], pmax[1])
            ylims!(ax, pmin[2], pmax[2])
            zlims!(ax, pmin[3], pmax[3])
            lines!(ax, gvecs)
            arrows3d!(ax, [Point3f(0, 0, 0)], @lift([$gvec]))
        end
    end

    # Note: `uobs[]` must already be updated before `tobs`
    _attenuation = Point2f[]
    attenuation = @lift push!(_attenuation, Point2f($tobs, abs(sum(M * uobs[])) / S0))

    ax = Axis(fig[2, 1]; xlabel = "t [µs]", title = "Signal attenuation")
    xlims!(ax, 0, TE)
    ylims!(ax, 0, 1.1)
    lines!(ax, attenuation)

    # Magnetization: magnitude and phase
    for (row, part, title) in ((1, :abs, "Magnitude"), (2, :phase, "Phase"))
        ax = mesh_axis(fig[row, 2], dim; title)
        colorrange = part == :abs ? (zero(umax), umax) : nothing
        plot = plot_field!(ax, femesh, uobs; part, colorrange)
        Colorbar(fig[row, 3], plot; label = field_label(part))
    end
    colsize!(fig.layout, 2, Relative(0.5))

    p.t = tobs
    p.u = uobs
    p.fig = fig

    # The figure is shown in a window with interactive backends only (other backends, such
    # as CairoMakie, would open a static image in a viewer)
    backend = Makie.current_backend()
    backend isa Module && nameof(backend) in (:GLMakie, :WGLMakie) && display(fig)

    p
end

function update_plotter!(p::Plotter, problem, gradient, u, t)
    if p.n % p.nupdate == 0
        # u updates first: important!
        p.u[] = u
        p.t[] = t
        sleep(0.01)
    end
    p.n += 1
end
