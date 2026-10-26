# Checks of computed values against exact identities and analytical solutions

γ = 2.67513e-4

@testset "Permeable plates vs analytical solution" begin
    # For a gradient across the plates, the 2D plates and the 3D slabs reduce to
    # the one-dimensional problem of the analytical matrix formalism
    profile = PGSE(2000.0, 6000.0)
    for setup ∈ [
        PlateSetup(; depth = 20.0, widths = fill(5.0, 4), refinement = 0.5),
        SlabSetup(; depth = 20.0, widths = fill(5.0, 4), height = 5.0, refinement = 2.0),
    ]
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
        matrices = assemble_matrices(model)
        dim = size(mesh.points[1], 1)
        dir = [1.0; zeros(dim - 1)]
        gradient = ScalarGradient(dir, profile, √(1000 / int_F²(profile)) / γ)
        ξ = solve(
            BTPDE(; model, matrices),
            gradient,
            IntervalConstantSolver(; timestep = 5.0),
        )
        volumes = get_cmpt_volumes(mesh)
        signal = compute_signal(matrices.M, ξ) / sum(volumes)
        acoeffs = analytical_coefficients(setup, coeffs)
        laplace = AnalyticalLaplace(; acoeffs..., eiglim = length2eig(0.3, 0.002), eigstep = 1e-8)
        amf = AnalyticalMatrixFormalism(;
            analytical_laplace = laplace,
            lap_mat = solve(laplace),
            volumes,
        )
        signal_analytical = solve(amf, gradient) / sum(volumes)
        @test real(signal) ≈ real(signal_analytical) rtol = 1e-2
    end
end
