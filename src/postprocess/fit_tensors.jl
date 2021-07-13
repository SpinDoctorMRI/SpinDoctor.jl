"""
    fit_tensors(directions, adc)

Fit effective diffusion tensors to directionalized ADCs.
The six components of the symmetric diffusion tensors are fitted by least
squares to the gradient directions and resulting ADCs.
"""
function fit_tensors(directions, adc)

    # Sizes
    ncompartment, nsequence, ndirection = size(adc)

    # Get vectors of gradient interactions
    g = directions
    G = zeros(ndirection, 6)
    G[:, 1] = g[1, :] .^ 2
    G[:, 2] = g[2, :] .^ 2
    G[:, 3] = g[3, :] .^ 2
    G[:, 4] = 2 .* g[1, :] .* g[2, :]
    G[:, 5] = 2 .* g[1, :] .* g[3, :]
    G[:, 6] = 2 .* g[2, :] .* g[3, :]

    tensors = Array{SMatrix{3,3},2}(undef, ncompartment, nsequence)
    for icmpt = 1:ncompartment, iseq = 1:nsequence

        # Deduce vector of effective diffusion tensor components usinag a least squares fit
        # of the directions to the computed ADCs
        Dvec = G \ adc[icmpt, iseq, :]

        # Tensorize component vectors (symmetric)
        tensors[icmpt, iseq] = SMatrix{3,3}(
            [
                Dvec[1, :] Dvec[4, :] Dvec[5, :]
                Dvec[4, :] Dvec[2, :] Dvec[6, :]
                Dvec[5, :] Dvec[6, :] Dvec[3, :]
            ],
        )

    end

    tensors
end
