# Dense random packings of disks and spheres in periodic boxes

"""
    periodic_packing(
        radii;
        dim = 2,
        fraction,
        gap = 0.0,
        rng = Random.default_rng(),
        maxiter = 200_000,
    )

Pack disks (`dim = 2`) or spheres (`dim = 3`) with the given `radii` in a periodic box,
without overlaps and with gaps of at least `gap` between neighbors. The box size is chosen
so that the cells fill the volume fraction `fraction` of the box.

The positions are found by a collective rearrangement: starting from random positions, the
radii grow towards their final values while the overlaps are removed by moving overlapping
pairs apart. Random packings of polydisperse disks jam at fractions of about 0.84 and
spheres at about 0.64 (counting the gaps as part of the cells); higher fractions cannot be
reached. If the packing jams (or takes more than `maxiter` relaxation sweeps) within one
percent of the final radii, the radii are reduced accordingly (with a warning); otherwise
an error is thrown.

Returns a cell configuration `(; centers, radii, box)`, with `centers` of size
`(dim, ncell)` inside `box = (lower, upper)`, for [`create_geometry`](@ref) with a
[`PeriodicBoxECS`](@ref).
"""
function periodic_packing(
    radii;
    dim = 2,
    fraction,
    gap = 0.0,
    rng = Random.default_rng(),
    maxiter = 200_000,
)
    0 < fraction < 1 || error("The fraction must be between 0 and 1")
    radii = collect(float.(radii))
    n = length(radii)
    T = eltype(radii)
    cellvolume(r) = dim == 2 ? pi * r^2 : 4 / 3 * pi * r^3
    L = (sum(cellvolume, radii) / fraction)^(1 / dim)

    # Exclusion radii (with half the gap on each side)
    R = radii .+ gap / 2
    2 * maximum(R) < L || error("The largest cell does not fit in the box")
    x = rand(rng, T, dim, n) .* L
    dx = zeros(T, dim, n)

    # Relax the overlaps at the radii s * R, then grow s towards one
    s = T(0.5)
    ds = T(0.02)
    xbest = copy(x)
    sbest = zero(T)
    # Overlaps of the exclusion radii smaller than a tenth of the gap leave a gap of at least
    # 0.9 gap between the cells
    tol = gap > 0 ? T(gap / 10) : T(1e-6) * minimum(R)
    iter = 0
    shakes = 0
    while true
        relaxed = false
        for _ = 1:1000
            iter += 1
            iter > maxiter && break
            overlap = remove_overlaps!(x, dx, s * R, L, tol)
            if overlap < tol
                relaxed = true
                break
            end
        end
        if relaxed
            xbest .= x
            sbest = s
            s == 1 && break
            s = min(one(T), s * (1 + ds))
            shakes = 0
        elseif shakes < 10 && iter <= maxiter
            # Stuck: shake the configuration
            x .+= (gap > 0 ? gap : minimum(R) / 10) / 2 .* (rand(rng, T, dim, n) .- T(0.5))
            shakes += 1
        else
            # Jammed at this size: grow more slowly from the last packing without overlaps
            ds /= 2
            shakes = 0
            x .= xbest
            s = min(one(T), sbest * (1 + ds))
            if ds < 1e-6 || iter > maxiter
                # Keep the last packing without overlaps, with smaller cells
                x .= xbest
                radii = sbest .* R .- gap / 2
                fraction_reached = sum(cellvolume, radii) / L^dim
                if sbest < 0.99
                    error("The packing jammed at a fraction of $fraction_reached")
                end
                @warn "The packing jammed: fraction $fraction_reached instead of $fraction"
                break
            end
        end
    end

    (;
        centers = mod.(x, L),
        radii = reshape(radii, 1, :),
        box = (zeros(T, dim), fill(L, dim)),
    )
end

# Move all overlapping pairs apart along their center lines, unless all overlaps are below
# `tol`, and return the largest overlap before moving. The pairs are the nearest periodic
# images in the box of size L, and also the other neighboring images for cells larger than
# a quarter of the box.
function remove_overlaps!(x, dx, R, L, tol)
    dim, n = size(x)
    T = eltype(x)
    fill!(dx, 0)
    maxoverlap = zero(T)
    shifts = vec(collect(Iterators.product(ntuple(_ -> -1:1, dim)...)))
    for j = 2:n, i = 1:(j-1)
        v0 = SVector{dim,T}(ntuple(k -> x[k, j] - x[k, i], dim))
        v0 = v0 .- L .* round.(v0 ./ L)
        for s in (R[i] + R[j] > L / 2 ? shifts : (ntuple(_ -> 0, dim),))
            v = v0 .+ L .* SVector{dim,T}(s)
            d = norm(v)
            overlap = R[i] + R[j] - d
            overlap > 0 || continue
            maxoverlap = max(maxoverlap, overlap)
            u = d > 0 ? v / d : SVector{dim,T}(ntuple(k -> k == 1, dim))
            dx[:, i] .-= overlap / 2 .* u
            dx[:, j] .+= overlap / 2 .* u
        end
    end
    if maxoverlap >= tol
        x .+= dx
        x .= mod.(x, L)
    end
    maxoverlap
end
