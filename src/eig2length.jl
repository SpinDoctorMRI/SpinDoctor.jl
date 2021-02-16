"""
    eig2length(λ, σ)

Convert Laplace eigenvalue to length scale.
"""
eig2length(λ, σ) = λ > 0 ? π * √(σ / λ) : Inf

"""
    length2eig(length, σ)

Convert length scale to Laplace eigenvalue.
"""
length2eig(length, σ) = σ * (π / length)^2
