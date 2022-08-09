function compute_JY(z, n, dim)
    if dim == 1
        J = cos.(z)
        Y = sin.(z)
        dJ = -sin.(z)
        dY = cos.(z)
    elseif dim == 2
        J = besselj.(n, z)
        Y = bessely.(n, z)
        dJ = (besselj.(n - 1, z) - besselj.(n + 1, z)) / 2
        dY = (bessely.(n - 1, z) - bessely.(n + 1, z)) / 2
    elseif dim == 3
        J = sphericalbesselj.(n, z)
        Y = sphericalbessely.(n, z)

        dJ = [copy(z);]
        dY = [copy(z);]
        ind = z .> 0
        zz = z[ind]
        ind = [ind;]
        @. dJ[ind] =
            (
                sphericalbesselj(n - 1, zz) - sphericalbesselj(n + 1, zz) -
                sphericalbesselj(n, zz) / zz
            ) / 2
        @. dY[ind] =
            (
                sphericalbessely(n - 1, zz) - sphericalbessely(n + 1, zz) -
                sphericalbessely(n, zz) / zz
            ) / 2

        ind = [z .== 0;]
        dJ[ind] .= (n == 1) / 3
        dY[ind] .= Inf

        if size(z) == ()
            J, Y, dJ, dY = J[1], Y[1], dJ[1], dY[1]
        end
    end

    J, Y, dJ, dY
end
