"""
    limit_lengthscale(lap_eig, λ_max)

Only keep modes with length scales larger than minimum
"""
function limit_lengthscale(lap_eig, λ_max)
    (; values, funcs, moments, massrelax) = lap_eig
    inds = values .≤ λ_max
    values = values[inds]
    funcs = funcs[:, inds]
    moments = [m[inds, inds] for m ∈ moments]
    massrelax = massrelax[inds, inds]
    (; values, funcs, moments, massrelax)
end
