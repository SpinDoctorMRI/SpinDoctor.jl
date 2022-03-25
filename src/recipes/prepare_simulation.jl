"""
    prepare_simulation(setup::AbstractSetup{T}; recreate::Bool = true)

"""
function prepare_simulation(setup::AbstractSetup{T}; recreate::Bool=false) where {T}
    # Get compartimentalized coefficient vectors
    coeffs = coefficients(setup)
    (; ρ, D, T₂, κ, γ) = coeffs

    mesh, surfaces, cells = @time create_geometry(setup; recreate=recreate)

    # We may also compute some useful quantities, including a scalar diffusion coefficient from
    # the diffusion tensors.
    volumes = get_cmpt_volumes(mesh)
    D_avg = 1 / 3 * tr.(D)' * volumes / sum(volumes)
    ncompartment = length(mesh.points)
    @info "Number of nodes per compartment:" size.(mesh.points, 2) # because legnth = npts * ndim

    model = Model{T}(; mesh, ρ, D, T₂, κ, γ, D_avg, volumes, ncompartment)

    ## Assemble finite element matrices
    matrices = @time assemble_matrices(model)

    model, matrices, surfaces, cells
end