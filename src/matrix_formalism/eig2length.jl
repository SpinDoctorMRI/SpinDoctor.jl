"""
    eig2length(λ, D)

Convert Laplace eigenvalue to length scale.
"""
eig2length(λ, D) = λ > 0 ? π * √(D / λ) : Inf

"""
    length2eig(length, D)

Convert length scale to Laplace eigenvalue.
"""
length2eig(length, D) = D * (π / length)^2
