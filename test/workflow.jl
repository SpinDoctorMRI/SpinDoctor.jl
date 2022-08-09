# TODO: Make tests pass for
# - T ∈ [Float32, Float64]
# - Setup ∈ [PlateSetup, DiskSetup, ExtrusionSetup, Setup, SphereSetup,
#            NeuronSetup]
# - Gradient ∈ [GeneralGradient, ScalarGradient] (all time profiles)

# Floating point type for simulations
T = Float64

@testset "Setup recipes" begin
    @testset "PlateSetup" begin
        setup = get_setup(PlateSetup{T})
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "SlabSetup" begin
        setup = get_setup(SlabSetup{T})
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "DiskSetup" begin
        setup = get_setup(DiskSetup{T})
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "CylinderSetup" begin
        setup = get_setup(CylinderSetup{T})
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "SphereSetup" begin
        setup = get_setup(SphereSetup{T})
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
    end

    @testset "NeuronSetup" begin
        setup = get_setup(NeuronSetup{T})
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup; meshdir = "meshfiles/neuron")
        model = Model(; mesh, coeffs...)
    end

end

# Geometrical setup
setup = get_setup(SlabSetup{T})
coeffs = get_coeffs(setup)

# Get compartimentalized coefficient vectors
mesh, = create_geometry(setup; recreate = true)
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
g⃗(t) = 1.0 * R * [sin(2π * t / TE), sin(20π * t / TE) / 5, cos(2π * t / TE)]
general_gradient = GeneralGradient{T,typeof(g⃗)}(; g⃗, TE)


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
    hadc = HADC(; model, matrices)
    @test_throws MethodError solve(hadc, general_gradient)
    @test solve(hadc, pgse_gradient) isa Vector{T}
    @test solve(hadc, ogse_gradient) isa Vector{T}
end

@testset "BTPDE" begin
    btpde = BTPDE(; model, matrices)
    @test solve(btpde, general_gradient) isa Vector{Complex{T}}
    @test solve(btpde, ogse_gradient) isa Vector{Complex{T}}
    solver = IntervalConstantSolver{T}(; θ = 0.5, timestep = 5)
    @test_throws MethodError solve(btpde, general_gradient, solver) isa Vector{Complex{T}}
    @test_throws ErrorException solve(btpde, ogse_gradient, solver) isa Vector{Complex{T}}
    @test solve(btpde, pgse_gradient, solver) isa Vector{Complex{T}}
end

@testset "Karger" begin
    directions = unitsphere(10)
    gradients = [
        ScalarGradient(collect(d), pgse_gradient.profile, pgse_gradient.amplitude) for
        d ∈ eachcol(directions)
    ]
    hadc = HADC(; model, matrices)
    adcs, = solve_multigrad(hadc, gradients)
    difftensors = fit_tensors(directions, adcs)
    karger = Karger(; model, difftensors)
    signal = solve(karger, pgse_gradient; timestep = 5.0)
end

@testset "Matrix formalism" begin
    laplace = Laplace(; model, matrices, neig_max = 400)
    lap_eig = solve(laplace)
    length_scales = eig2length.(lap_eig.values, D_avg)
    length_scale = 3
    λ_max = length2eig(length_scale, D_avg)
    lap_eig = limit_lengthscale(lap_eig, λ_max)
    @test all(λ -> eig2length(λ, D_avg) ≥ length_scale, lap_eig.values)
    mf = MatrixFormalism(; model, matrices, lap_eig)
    @test solve(mf, general_gradient) isa Vector{Complex{T}}
    @test solve(mf, pgse_gradient) isa Vector{Complex{T}}
    @test solve(mf, ogse_gradient) isa Vector{Complex{T}}
end

@testset "Analytical" begin
    length_scale = 1.0
    eigstep = 1e-8
    eiglim = length2eig(length_scale, D_avg)
    analytical_coeffs = analytical_coefficients(setup, coeffs)
    analytical_laplace = AnalyticalLaplace(; analytical_coeffs..., eiglim, eigstep)
    lap_mat = solve(analytical_laplace)
    analytical_mf = AnalyticalMatrixFormalism(; analytical_laplace, lap_mat, volumes)
    @test_throws MethodError solve(analytical_mf, general_gradient)
    @test_throws ErrorException solve(analytical_mf, ogse_gradient)
    @test solve(analytical_mf, pgse_gradient) isa Complex{T}
end
