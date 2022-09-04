function compute_β(params, α, n)
    I = compute_int_I(params, α, α, n)
    1 / √(2 * sum(I))
end
