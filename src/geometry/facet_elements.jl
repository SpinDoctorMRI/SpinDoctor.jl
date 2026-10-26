"""
    facet_elements(elements, facets)

For each facet (column of `facets`), find an element (column of `elements`) that has the
facet as one of its faces. The result is zero for facets that are not a face of any element.
"""
facet_elements(elements, facets) = facet_elements(elements, facets, Val(size(facets, 1)))

function facet_elements(elements, facets, ::Val{d}) where {d}
    nbasis = size(elements, 1)
    face(f, j) = sort(SVector{d}(ntuple(k -> f[k, j], Val(d))))

    # Look up the faces of all elements (drop one vertex at a time) among the facets
    lookup = Dict(face(facets, j) => j for j ∈ axes(facets, 2))
    ielements = zeros(Int, size(facets, 2))
    for e ∈ axes(elements, 2), skip = 1:nbasis
        f = sort(SVector{d}(ntuple(k -> elements[k < skip ? k : k + 1, e], Val(d))))
        j = get(lookup, f, 0)
        j > 0 && (ielements[j] = e)
    end

    # Repeated facets share the element of the one in the lookup table
    for j ∈ axes(facets, 2)
        ielements[j] = ielements[lookup[face(facets, j)]]
    end

    ielements
end
