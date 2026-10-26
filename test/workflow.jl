# TODO: Make tests pass for
# - T ∈ [Float32, Float64]
# - Setup ∈ [PlateSetup, DiskSetup, ExtrusionSetup, SlabSetup, CylinderSetup,
#       SphereSetup, NeuronSetup]
# - Gradient ∈ [GeneralGradient, ScalarGradient] (all time profiles)

T = Float64

@testset "Setup recipes" begin
    @testset "PlateSetup" begin
        setup = PlateSetup()
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "SlabSetup" begin
        setup = SlabSetup()
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "DiskSetup" begin
        setup = DiskSetup()
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "CylinderSetup" begin
        setup = CylinderSetup(; ecs = ConvexHullECS(; margin = 0.5))
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "SphereSetup" begin
        setup = SphereSetup()
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "NeuronSetup" begin
        setup = NeuronSetup()
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup; savedir = "meshfiles/neuron")
        model = Model(; mesh, coeffs...)
    end
end

# Geometrical setup
setup = SlabSetup(; depth = 20.0, widths = fill(5.0, 4), height = 20.0, refinement = 10.0)
coeffs = get_coeffs(setup)

# Get compartimentalized coefficient vectors
mesh, surfaces, cells = create_geometry(setup)
model = Model(; mesh, coeffs...)
volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)
ncompartment = length(model.mesh.points)
matrices = assemble_matrices(model);

## Gradients
dir = [1.0, 1.0, 1.0]
profile = PGSE(2000.0, 6000.0)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
pgse_gradient = ScalarGradient(dir, profile, g)

dir = [1.0, 1.0, 1.0]
profile = CosOGSE(5000.0, 5000.0, 2)
b = 1000
g = √(b / int_F²(profile)) / coeffs.γ
ogse_gradient = ScalarGradient(dir, profile, g)

TE = 5000
φ = -π / 6
R = [cos(φ) sin(φ) 0; -sin(φ) cos(φ) 0; 0 0 1]
gvec(t) = 1.0 * R * [sin(2π * t / TE), sin(20π * t / TE) / 5, cos(2π * t / TE)]
general_gradient = GeneralGradient(; gvec, TE)

@testset "Signal" begin
    ρ = initial_conditions(model)
    @test compute_signal(matrices.M, ρ) isa Complex{T}
    @test compute_signal.(matrices.M_cmpts, split_field(model.mesh, ρ)) isa
          Vector{Complex{T}}
end

@testset "ADC STA" begin
    @test compute_adc_sta(model, pgse_gradient) isa Vector{T}
    @test compute_adc_sta(model, ogse_gradient) isa Vector{T}
    @test_throws MethodError compute_adc_sta(model, general_gradient)
end

@testset "HADC" begin
    @test_throws MethodError solve_hadc(model, matrices, general_gradient)
    @test solve_hadc(model, matrices, pgse_gradient) isa Vector{T}
    @test solve_hadc(model, matrices, ogse_gradient) isa Vector{T}
end

@testset "BTPDE" begin
    @test solve_btpde(model, matrices, general_gradient) isa Vector{Complex{T}}
    @test solve_btpde(model, matrices, ogse_gradient) isa Vector{Complex{T}}
    solver = IntervalConstantSolver(; θ = 0.5, timestep = 5.0)
    @test_throws MethodError solve_btpde(model, matrices, general_gradient, solver)
    @test_throws ErrorException solve_btpde(model, matrices, ogse_gradient, solver)
    @test solve_btpde(model, matrices, pgse_gradient, solver) isa Vector{Complex{T}}
end

@testset "Karger" begin
    directions = unitsphere(10)
    gradients = [
        ScalarGradient(collect(d), pgse_gradient.profile, pgse_gradient.amplitude) for
        d ∈ eachcol(directions)
    ]
    adcs = [solve_hadc(model, matrices, grad) for grad ∈ gradients]
    difftensors = fit_tensors(directions, adcs)
    signal = solve_karger(model, difftensors, pgse_gradient; timestep = 5.0)
    @test signal isa Vector{Complex{T}}
end

@testset "Matrix formalism" begin
    lap_eig = solve_laplace(model, matrices; neig_max = 400)
    length_scales = eig2length.(lap_eig.values, D_avg)
    length_scale = 3
    λ_max = length2eig(length_scale, D_avg)
    lap_eig = limit_lengthscale(lap_eig, λ_max)
    @test all(λ -> eig2length(λ, D_avg) ≥ length_scale, lap_eig.values)
    @test solve_mf(model, matrices, lap_eig, general_gradient) isa Vector{Complex{T}}
    @test solve_mf(model, matrices, lap_eig, pgse_gradient) isa Vector{Complex{T}}
    @test solve_mf(model, matrices, lap_eig, ogse_gradient) isa Vector{Complex{T}}
end

@testset "Analytical" begin
    length_scale = 1.0
    eigstep = 1e-8
    eiglim = length2eig(length_scale, D_avg)
    analytical_coeffs = analytical_coefficients(setup, coeffs)
    laplace = solve_analytical_laplace(; analytical_coeffs..., eiglim, eigstep)
    @test_throws MethodError solve_analytical_mf(laplace, volumes, general_gradient)
    @test_throws ErrorException solve_analytical_mf(laplace, volumes, ogse_gradient)
    @test solve_analytical_mf(laplace, volumes, pgse_gradient) isa Complex{T}
end
