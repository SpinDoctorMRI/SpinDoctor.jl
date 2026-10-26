"""
    facet_elements(elements, facets)

For each facet (column of `facets`), find an element (column of `elements`) that has the
facet as one of its faces. The result is zero for facets that are not a face of any element.
"""
facet_elements(elements, facets) = facet_elements(elements, facets, Val(size(facets, 1)))

function facet_elements(elements, facets, ::Val{d}) where {d}
    nbasis = size(elements, 1)
    face(f, j) = sort(SVector{d}(ntuple(k -> f[k, j], d)))

    # All faces of all elements: drop one vertex at a time
    lookup = Dict{SVector{d,Int},Int}()
    sizehint!(lookup, nbasis * size(elements, 2))
    for e ∈ axes(elements, 2), skip = 1:nbasis
        f = sort(SVector{d}(ntuple(k -> elements[k < skip ? k : k + 1, e], d)))
        lookup[f] = e
    end

    [get(lookup, face(facets, j), 0) for j ∈ axes(facets, 2)]
end
