"""
    assemble_flux_matrices(points, facets)

Assemble flux matrix (`Q`) for each compartment and boundary.
"""
function assemble_flux_matrices(points, facets)

    # Sizes
    ncompartment, nboundary = size(facets)

    flux_matrices = [
        assemble_flux_matrix(facets[icmpt, iboundary]', points[icmpt]') for
        icmpt = 1:ncompartment, iboundary = 1:nboundary
    ]

    flux_matrices
end
