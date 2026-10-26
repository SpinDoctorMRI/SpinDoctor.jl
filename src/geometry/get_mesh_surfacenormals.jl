"""
    get_mesh_surfacenormals(points, elements, facets)

Compute areas, centers and outward unit normals of the boundary `facets` of a mesh
compartment with `points` and `elements`.
"""
function get_mesh_surfacenormals(points, elements, facets)
    dim = size(points, 1)
    areas, facet_centers, normals = get_mesh_surface(points, facets)

    # The element of each facet lies on the inner side of the facet
    ielements = facet_elements(elements, facets)
    all(>(0), ielements) || error("Some facets are not faces of the elements")
    element_centers = reshape(mean(points[:, elements[:, ielements]]; dims = 2), dim, :)

    # Orient the normals from the element center towards the facet center
    orientations = sign.(sum((facet_centers - element_centers) .* normals; dims = 1))
    normals = orientations .* normals

    areas, facet_centers, normals
end
