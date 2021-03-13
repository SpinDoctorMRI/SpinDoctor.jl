"""
    get_cmpt_volumes(mesh)

Get compartment volumes.
"""
function get_cmpt_volumes(mesh)
    map((p, e) -> sum(get_mesh_volumes(p, e)[1]), mesh.points, mesh.elements)
end


"""
    get_mesh_volumes(points, elements)

Get total volume and volumes of each tetrahedron of mesh.
"""
function get_mesh_volumes(points, elements)

    # Sizes
    nelement = size(elements, 2)

    # Elements
    x = reshape(points[:, elements], 3, 4, nelement)

    # Element centers
    centers = mean(x, dims = 2)[:, 1, :]

    # Element volumes
    areavectors = [cross(x[:, 2, i] - x[:, 4, i], x[:, 3, i] - x[:, 4, i]) for i = 1:nelement]

    volumes = [1 / 6 * abs(dot(x[:, 1, i] - x[:, 4, i], areavectors[i])) for i = 1:nelement]

    volumes, centers
end
