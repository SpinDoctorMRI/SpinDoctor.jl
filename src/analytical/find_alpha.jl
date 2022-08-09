function find_α(params, α_min, α_max, dα)
    (; dim) = params

    if dim == 1
        α = find_zeros(α -> α_func(params, α, 0), α_min, α_max)
        n = zeros(Int, size(α))
    else
        α = Float64[]
        n = Int[]
        α₀ = [0.0]
        j = 0
        # pl = nothing
        # αα = LinRange(α_min, α_max, 1000)
        while !isempty(α₀) && α₀[1] < α_max
            @info "Computing zeros of Fₙ(α) for n = $j"
            α₀ = find_zeros(α -> α_func(params, α, j), α_min, α_max)
            append!(α, α₀)
            append!(n, fill(j, size(α₀)))
            # if j == 0
            #     pl = lines(αα,  α -> α_func(params, α, j))
            #     display(pl)
            # else
            #     lines!(αα, α -> α_func(params, α, j))
            # end
            j = j + 1
        end

        inds = sortperm(α)
        α = α[inds]
        n = n[inds]
    end

    # Number of layers
    L = length(params.D)

    inds_keep = dα / 100 .< α .< α_max
    α = α[inds_keep]
    n = n[inds_keep]

    if params.W[L] < 1e-12 && (length(params.W) == L || params.W[end] < 1e-12)
        # Add ground eigenmode explicitly
        pushfirst!(α, 0)
        pushfirst!(n, 0)
    end

    α, n
end
