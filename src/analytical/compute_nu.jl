function compute_ν(n, dim)
    if dim == 1
        ν = 0
    elseif dim == 2
        ν = n^2
    elseif dim == 3
        ν = n * (n + 1)
    end
    ν
end
