"""
    get_mesh_surfacenormals(points, elements, facets)

Computing outwards oriented surface normals on a connected mesh compartment.
"""
function get_mesh_surfacenormals(points, elements, facets)
    dim = size(points, 1)
    nfacet = size(facets, 2)
    nelement = size(elements, 2)

    _, areas, facet_centers, normals = get_mesh_surface(points, facets)

    # Find elements to which the facets belong
    if dim == 2
        combs = [[1, 2], [1, 3], [2, 3]]
    elseif dim == 3
        combs = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
    end
    ef = [SVector{dim}(sort(elements[comb, i])) for comb ∈ combs for i = 1:nelement]
    sort!(ef)
    ef_labels = repeat(1:nelement, dim + 1)

    f = SVector{dim}.(sort.(eachcol(facets)))
    matchinds = fill(0, length(f))
    for i ∈ eachindex(f)
        j = findfirst(isequal(f[i]), ef)
        matchinds[i] = ef_labels[j]
    end

    elements_match = elements[:, matchinds]

    # Identify centers of matching elements
    x = points[:, elements_match]
    element_centers = reshape(mean(x; dims = 2), dim, nfacet)

    # Create outward directed vectors (from cell center to outer facet)
    outwards_vectors = facet_centers - element_centers

    # Determine orientation of normals (+ for out, - for in)
    orientations = [sign(outwards_vectors[:, i]' * normals[:, i]) for i = 1:nfacet]

    # Orient all normals outwards
    normals = orientations' .* normals

    areas, facet_centers, normals
end
