function compute_ν(n, dim)
    if dim == 2
        ν = n^2
    else
        ν = n * (n + 1)
    end

    ν
end