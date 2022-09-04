function α_func(params, α, n)
    (; r, D, W, dim) = params

    m = length(D)
    bc1 = zeros(size(α))
    bc2 = zeros(size(α))

    for i = 1:length(α)
        bc = compute_bc(params, α[i], n)
        bc1[i] = bc[m, 1]
        bc2[i] = bc[m, 2]
    end

    rD = r[m] / √D[m]

    if dim == 1
        f1 = -sqrt(D[m]) * α .* sin(α * rD) + W[m] * cos(α * rD)
        f2 = sqrt(D[m]) * α .* cos(α * rD) + W[m] * sin(α * rD)
    elseif dim == 2
        f1 =
            √D[m] * α .* (besselj(n - 1, α * rD) - besselj(n + 1, α * rD)) / 2 +
            W[m] * besselj(n, α * rD)
        f2 =
            √D[m] * α .* (bessely(n - 1, α * rD) - bessely(n + 1, α * rD)) / 2 +
            W[m] * bessely(n, α * rD)
    elseif dim == 3
        f1 =
            √D[m] * α .* (
                sphericalbesselj(n - 1, α * rD) - sphericalbesselj(n + 1, α * rD) -
                sphericalbesselj(n, α * rD) ./ (α * rD)
            ) / 2 + W[m] * sphericalbesselj(n, α * rD)
        f2 =
            √D[m] * α .* (
                sphericalbessely(n - 1, α * rD) - sphericalbessely(n + 1, α * rD) -
                sphericalbessely(n, α * rD) ./ (α * rD)
            ) / 2 + W[m] * sphericalbessely(n, α * rD)
    end

    bc1 .* f1 + bc2 .* f2
end
