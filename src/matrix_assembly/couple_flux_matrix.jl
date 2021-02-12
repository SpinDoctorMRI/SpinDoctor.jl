# Create global flux matrix with coupling between compartments
function couple_flux_matrix(mesh, Q)

    # Extract mesh fields
    @unpack ncmpt, nboundary, cmpt_inds, points, facets, elements, boundary_markers = mesh

    # Sizes (not the same as `mesh.cmpt_inds`)
    inds_cmpts = cumsum([0; size.(points, 2)])
    get_inds(icmpt) = inds_cmpts[icmpt]+1:inds_cmpts[icmpt+1]

    # Assemble sparse block diagonal flux matrix
    couple_Q = blockdiag(Q...)

    # Couple flux matrix at boundaries
    for iboundary = 1:nboundary

        # Check if boundary is an interface between two compartments
        cmpts_touch = findall(boundary_markers[:, iboundary])
        if length(cmpts_touch) == 2
            # Extract flux matrices from corresponding compartments
            cmpt1, cmpt2 = cmpts_touch[1], cmpts_touch[2]
            Q1 = Q[cmpt1]
            Q2 = Q[cmpt2]
            Q12 = spzeros(size(Q1, 1), size(Q2, 2))
            Q21 = spzeros(size(Q2, 1), size(Q1, 2))

            # Identify pairs of points using global numbering
            inds1 = unique(facets[cmpt1, iboundary])
            inds2 = unique(facets[cmpt2, iboundary])
            nind = length(inds1)
            pairs = cmpt_inds[cmpt1][inds1] .== cmpt_inds[cmpt2][inds2]'

            # Couple matrices
            for i = 1:nind
                # Loacal incdices of coupled points
                i1 = inds1[i]
                i2 = inds2[findfirst(pairs[i, :])]

                # Create flux coupling between points
                Q12[:, i2] = -Q1[:, i1]
                Q21[:, i1] = -Q2[:, i2]
            end

            # Store compartment interface coupling in global matrix
            couple_Q[get_inds(cmpt1), get_inds(cmpt2)] = Q12
            couple_Q[get_inds(cmpt2), get_inds(cmpt1)] = Q21
        end
    end

    couple_Q
end
