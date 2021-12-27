"""
    fit_tensors(directions, adcs)

Fit effective diffusion tensors to directionalized ADCs.
The six components of the symmetric diffusion tensors are fitted by least
squares to the gradient directions and resulting ADCs.
"""
function fit_tensors(directions, adcs)
    T = eltype(directions)

    adc = reduce(hcat, adcs)
    ncompartment, ndirection = size(adc)

    # Get vectors of gradient interactions
    d = directions
    G = zeros(ndirection, 6)
    G[:, 1] = d[1, :] .^ 2
    G[:, 2] = d[2, :] .^ 2
    G[:, 3] = d[3, :] .^ 2
    G[:, 4] = 2 .* d[1, :] .* d[2, :]
    G[:, 5] = 2 .* d[1, :] .* d[3, :]
    G[:, 6] = 2 .* d[2, :] .* d[3, :]

    tensors = Vector{SMatrix{3,3,T,9}}(undef, ncompartment)
    for icmpt = 1:ncompartment
        # Deduce vector of effective diffusion tensor components usinag a least squares fit
        # of the directions to the computed ADCs
        Dvec = G \ adc[icmpt, :]

        # Tensorize component vectors (symmetric)
        tensors[icmpt] = SMatrix{3,3}(
            [
                Dvec[1, :] Dvec[4, :] Dvec[5, :]
                Dvec[4, :] Dvec[2, :] Dvec[6, :]
                Dvec[5, :] Dvec[6, :] Dvec[3, :]
            ],
        )

    end

    tensors
end
