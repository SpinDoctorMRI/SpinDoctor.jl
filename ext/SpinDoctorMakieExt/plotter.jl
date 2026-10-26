function initialize_plotter!(p::Plotter, problem, gradient, u, t)
    T = real(eltype(u))
    p.n = 1
    tobs = Observable(t)
    uobs = Observable(u)

    femesh = problem.model.mesh
    M = problem.matrices.M
    ncompartment, nboundary = size(femesh.facets)
    dim = size(femesh.points[1], 1)

    umax = maximum(abs, u)
    S0 = abs(sum(M * u))
    TE = echotime(gradient)

    npoint_cmpts = size.(femesh.points, 2)
    inds = [0; cumsum(npoint_cmpts[1:end])]
    u_cmpts = [@lift($uobs[(1+inds[i]):inds[i+1]]) for i = 1:ncompartment]

    fig = Figure()

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

    gm = fig[1, 2] = GridLayout()
    colorrange = (0, umax)
    if dim == 2
        ax = Axis(gm[1, 1]; title = "Magnetization (magnitude)")
        for icmpt = 1:ncompartment
            elements = femesh.elements[icmpt]
            points = femesh.points[icmpt]
            color = @lift abs.($(u_cmpts[icmpt]))
            mesh!(ax, points', elements'; color, shading = NoShading, colorrange)
        end
    elseif dim == 3
        ax = Axis3(gm[1, 1]; title = "Magnetization (magnitude)")
        ax.aspect = :data
        for icmpt = 1:ncompartment
            color = @lift abs.($(u_cmpts[icmpt]))
            points = femesh.points[icmpt]
            for iboundary = 1:nboundary
                facets = femesh.facets[icmpt, iboundary]
                if !isempty(facets)
                    mesh!(ax, points', facets'; color, shading = NoShading, colorrange)
                end
            end
        end
    end
    Colorbar(gm[1, 2]; limits = colorrange)

    gm = fig[2, 2] = GridLayout()
    colormap = :cyclic_mrybm_35_75_c68_n256
    colorrange = (-pi, pi)
    if dim == 2
        ax = Axis(gm[1, 1]; title = "Magnetization (phase-shift)")
        for icmpt = 1:ncompartment
            color = @lift angle.($(u_cmpts[icmpt]))
            elements = femesh.elements[icmpt]
            points = femesh.points[icmpt]
            mesh!(ax, points', elements'; color, shading = NoShading, colorrange, colormap)
        end
    elseif dim == 3
        ax = Axis3(gm[1, 1]; title = "Magnetization (phase-shift)")
        ax.aspect = :data
        for icmpt = 1:ncompartment
            color = @lift angle.($(u_cmpts[icmpt]))
            points = femesh.points[icmpt]
            for iboundary = 1:nboundary
                facets = femesh.facets[icmpt, iboundary]
                if !isempty(facets)
                    mesh!(
                        ax,
                        points',
                        facets';
                        color,
                        shading = NoShading,
                        colorrange,
                        colormap,
                    )
                end
            end
        end
    end
    Colorbar(gm[1, 2]; limits = colorrange, colormap)

    p.t = tobs
    p.u = uobs
    p.fig = fig
    display(fig)

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
