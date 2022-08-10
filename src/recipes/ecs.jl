abstract type AbstractECS end

struct NoECS <: AbstractECS end

@with_kw struct BoxECS{T} <: AbstractECS
    margin::T = 1.0
    @assert margin > 0
end

@with_kw struct ConvexHullECS{T} <: AbstractECS
    margin::T = 1.0
    @assert margin > 0
    threshold::T = π / 20
    @assert threshold ≥ 0
    maxiter::Int = 10
    @assert maxiter ≥ 0
end

@with_kw struct TightWrapECS{T} <: AbstractECS
    margin::T = 1.0
    @assert margin > 0
end

function create_surfaces(::NoECS, points)
    T = eltype(points)
    dim = size(points, 1)
    (; points = zeros(T, dim, 0), facets = zeros(Int, dim, 0), region = zeros(T, dim, 0))
end

function create_surfaces(ecs::BoxECS{T}, points) where {T}
    dim = size(points, 1)
    (; margin) = ecs

    # Determine bounds of domain
    pmin = minimum(points; dims = 2) .- margin
    pmax = maximum(points; dims = 2) .+ margin

    if dim == 2
        # Define four corners of box ECS and their corresponding edges
        points = [
            pmin[1] pmax[1] pmax[1] pmin[1]
            pmin[2] pmin[2] pmax[2] pmax[2]
        ]
        facets = [
            1 2 3 4
            2 3 4 1
        ]
    elseif dim == 3
        # Define eight corners of box ECS and their corresponding edges
        points = [
            pmin[1] pmax[1] pmax[1] pmin[1] pmin[1] pmax[1] pmax[1] pmin[1]
            pmin[2] pmin[2] pmax[2] pmax[2] pmin[2] pmin[2] pmax[2] pmax[2]
            pmin[3] pmin[3] pmin[3] pmin[3] pmax[3] pmax[3] pmax[3] pmax[3]
        ]
        facets = [
            1 1 1 1 2 2 3 3 4 4 5 5
            2 3 2 6 3 7 4 8 1 5 6 7
            3 4 6 5 7 6 8 7 5 8 7 8
        ]
    end

    # The point half a margin in from a corner is guaranteed to be inside ECS
    region = pmin .+ margin / 2

    (; points, facets, region)
end

function create_surfaces(ecs::ConvexHullECS, points)
    (; margin, threshold, maxiter) = ecs
    threshold_inv = 1 / threshold
    dim = size(points, 1)

    facets, points = convexhull(points)
    nfacet = size(facets, 2)

    # Average point, is inside the convex hull by definition (convex
    # combination of points)
    c = sum(points; dims = 2) / size(points, 2)

    # Get outward pointing normals
    _, facet_centers, normals = get_mesh_surface(points, facets)
    o = sum((c .- facet_centers) .* normals; dims = 1)
    orientations = 1 .- 2 .* (o .> 0)
    normals .*= orientations

    # Add a point in the middle of each degenerate facet, and move it slightly
    # outwards
    if dim == 3
        for iter = 1:maxiter
            # Check for degenerate facets
            iadd = fill(false, nfacet)
            for i = 1:nfacet
                f = facets[:, i]
                a = points[:, f[2]] - points[:, f[1]]
                b = points[:, f[3]] - points[:, f[1]]
                c = points[:, f[3]] - points[:, f[2]]
                α = acos(dot(a, b) / norm(a) / norm(b))
                β = acos(dot(a, c) / norm(a) / norm(c))
                if abs(α) < threshold || abs(β) < threshold
                    iadd[i] = true
                end
            end

            padd = facet_centers[:, iadd] .+ 0.05 .* margin .* normals[:, iadd]
            points = [points padd]

            facets, points = convexhull(points)
            nfacet = size(facets, 2)

            # Average point, is inside the convex hull by definition (convex
            # combination of points)
            c = sum(points; dims = 2) / size(points, 2)

            # Get outward pointing normals
            _, facet_centers, normals = get_mesh_surface(points, facets)
            o = sum((c .- facet_centers) .* normals; dims = 1)
            orientations = 1 .- 2 .* (o .> 0)
            normals .*= orientations

            @info "Breaking up degenerate ECS facets" iter sum(iadd)
            sum(iadd) > 0 || break
        end
    end

    # Move each facet one margin outwards
    dirs = zero(points)
    for i = 1:size(facets, 2)
        # dirs[:, facets[:, i]] .+= normals[:, i]
        @views dirs[:, facets[:, i]] .+= normals[:, i]
    end
    ihull = map(c -> any(!≈(0), c), eachcol(dirs))
    points[:, ihull] .+=
        margin .* dirs[:, ihull] ./ sqrt.(sum(abs2, dirs[:, ihull]; dims = 1))

    # The point half a margin out from the first point is guaranteed inside
    # the ECS
    i = findfirst(ihull)
    region = points[:, i] - margin / 2 * dirs[:, i] / norm(dirs[:, i])

    (; points, facets, region)
end

function create_surfaces(ecs::TightWrapECS, points)
    error("Not implemented")
end
