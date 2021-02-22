""" Assemble flux matrix (Q) in compartment. """
function assemble_flux_matrix_cmpt(points, facets, κ)

    # Sizes
    npoint = size(points, 2)
    nboundary = length(facets)

    # Initialize output matrix
    flux_matrix = spzeros(size(points, 2), size(points, 2))

    # Add block in matrix for each touching boundary
    for iboundary = 1:nboundary

        # Check that there is coupling
        if κ[iboundary] > 1e-16

            # Extract boundary facets
            neumann = facets[iboundary]'

            # Only proceed if the boundary touches the compartment
            if !isempty(neumann)

                # Identify nodes in boundary
                neumann_nodes = unique(neumann)

                # Set weigths to boundary permeability for boundary nodes, else 0
                weights = zeros(npoint, 1)
                weights[neumann_nodes] .= κ[iboundary]

                # Add block to flux matrix (the blocks do not overlap)
                flux_matrix += assemble_flux_matrix(neumann, points', weights)
            end
        end
    end

    flux_matrix
end
