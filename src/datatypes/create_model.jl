function create_model(setup)
    coeffs = get_coefficients(setup)
    mesh, = create_geometry(setup)

    Model(; setup.name, mesh, coeffs...)
end
