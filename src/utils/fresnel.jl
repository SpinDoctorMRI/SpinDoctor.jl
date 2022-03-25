"""
    fresnels(x)

Compute the Fresnel S integral of type 0 (integral of `sin(πt²/2)`).
"""
fresnels(x) = (1 + im) / 4 * (erf((1 + im)√π / 2 * x) - im * erf((1 - im)√π / 2 * x))

"""
    fresnelc(x)

Compute the Fresnel C integral of type 0 (integral of `cos(πt²/2)`).
"""
fresnelc(x) = (1 - im) / 4 * (erf((1 + im)√π / 2 * x) + im * erf((1 - im)√π / 2 * x))
