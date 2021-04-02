"""
    fit_adc(b, S)

Fit Apparent Diffusion Coefficient (ADC) using a polynomial logfit of the normalized signal `S` against the b-values `b`.
"""
function fit_adc(b::Array{T,1}, S::Array{T,1}) where {T<:Real}

    # Fit log of signal
    p = fit(b, log.(S))

    # ADC is the negative linear term
    -p[1]
end
