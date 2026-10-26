# Tests of the GPU time stepping (CUDA extension), comparing with the CPU. They require an
# NVIDIA GPU, and are run separately from the main tests:
#
#     julia --project=test/gpu test/gpu/runtests.jl

using CUDA
using CUDSS
using LinearAlgebra
using Random
using SpinDoctor
using Test

CUDA.functional() || error("No functional CUDA GPU")
CUDA.allowscalar(false)

@testset "GPU" begin
    setup = CylinderSetup(;
        ncell = 3,
        rmin = 2.0,
        rmax = 4.0,
        layersizes = [0.6, 1.0],
        ecs = BoxECS(; margin = 2.0),
        height = 4.0,
    )
    nlayer = 2
    coeffs = coefficients(
        setup;
        D = (; cell = [0.002 * I(3), 0.0005 * I(3)], ecs = 0.002 * I(3)),
        T2 = (; cell = [Inf, Inf], ecs = Inf),
        rho = (; cell = [1.0, 1.0], ecs = 1.0),
        kappa = (;
            cell_interfaces = [1e-4],
            cell_boundaries = [0, 0],
            cell_ecs = 1e-4,
            ecs = 0,
        ),
        gamma = 2.67513e-4,
    )
    mesh, = create_geometry(setup; rng = Xoshiro(1))
    model = Model(; mesh, coeffs...)
    matrices = assemble_matrices(model)
    dir = normalize([1.0, 1.0, 0.5])

    @testset "BTPDE, constant gradient" begin
        gradient =
            ScalarGradient(dir, PGSE(2000.0, 6000.0); bvalue = 1000, gamma = model.gamma)
        cpu = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 50.0))
        printer = Printer(; nupdate = 100, verbosity = 0)
        gpu = solve_btpde(
            model,
            matrices,
            gradient,
            ThetaSolver(; timestep = 50.0, arraytype = CuArray);
            callbacks = [printer],
        )
        @test gpu isa Vector{ComplexF64}
        @test gpu ≈ cpu rtol = 1e-10
    end

    @testset "BTPDE, oscillating gradient" begin
        gradient = ScalarGradient(
            dir,
            CosOGSE(2000.0, 5000.0, 2);
            bvalue = 500,
            gamma = model.gamma,
        )
        cpu = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 100.0))
        gpu = solve_btpde(
            model,
            matrices,
            gradient,
            ThetaSolver(; timestep = 100.0, arraytype = CuArray),
        )
        @test gpu ≈ cpu rtol = 1e-8
    end

    @testset "HADC" begin
        profile = PGSE(2000.0, 6000.0)
        cpu = solve_hadc_tensors(model, matrices, profile; timestep = 50.0)
        gpu = solve_hadc_tensors(
            model,
            matrices,
            profile;
            timestep = 50.0,
            arraytype = CuArray,
        )
        @test all(gpu .≈ cpu)
    end

    @testset "Periodic boundary conditions" begin
        setup = CylinderSetup(;
            ncell = 8,
            rmin = 1.0,
            rmax = 2.0,
            nsidewall = 16,
            ecs = PeriodicBoxECS(; fraction = 0.6),
            height = 2.0,
        )
        coeffs = coefficients(
            setup;
            D = (; cell = [0.002 * I(3)], ecs = 0.002 * I(3)),
            T2 = (; cell = [Inf], ecs = Inf),
            rho = (; cell = [1.0], ecs = 1.0),
            kappa = (;
                cell_interfaces = Float64[],
                cell_boundaries = [0],
                cell_ecs = 1e-4,
                ecs = 0,
            ),
            gamma = 2.67513e-4,
        )
        mesh, = create_geometry(setup; rng = Xoshiro(1))
        model = Model(; mesh, coeffs...)
        matrices = assemble_matrices(model; periodic = true)
        profile = PGSE(2000.0, 6000.0)
        gradient = ScalarGradient(dir, profile; bvalue = 1000, gamma = model.gamma)
        cpu = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 50.0))
        solver = ThetaSolver(; timestep = 50.0, arraytype = CuArray)
        gpu = solve_btpde(model, matrices, gradient, solver)
        @test gpu ≈ cpu rtol = 1e-8
        cpu = solve_hadc_tensors(model, matrices, profile; timestep = 50.0)
        gpu = solve_hadc_tensors(
            model,
            matrices,
            profile;
            timestep = 50.0,
            arraytype = CuArray,
        )
        @test all(gpu .≈ cpu)
    end
end
