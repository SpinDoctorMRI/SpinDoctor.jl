"""
    fit_adc(bvalues, signals)

Fit the apparent diffusion coefficient (ADC) using a polynomial logfit of the normalized signals
`signals` against the b-values `bvalues`.
"""
function fit_adc(bvalues, signals)
    # Fit log of signal
    p = fit(bvalues, log.(abs.(signals)))

    # ADC is the negative linear term
    -p[1]
end
