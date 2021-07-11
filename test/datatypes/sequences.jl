@testset "sequences.jl" begin

    @testset "PGSE" begin
        δ, Δ = 500.0, 2000.0
        f = PGSE(δ, Δ)
        @test f.δ ≈ δ
        @test f.Δ ≈ Δ
        @test f(-100.0) ≈ 0.0
        @test f(prevfloat(0.0)) ≈ 0.0
        @test f(0.0) ≈ 1.0
        @test f(prevfloat(δ)) ≈ 1.0
        @test f(δ) ≈ 0.0
        @test f(nextfloat(δ)) ≈ 0.0
        @test f(prevfloat(Δ)) ≈ 0.0
        @test f(Δ) ≈ -1.0
        @test f(nextfloat(Δ)) ≈ -1.0
        @test f(Δ + δ + 100.0) ≈ 0.0

        # Integral
        T = LinRange(-δ, Δ + 2δ, 10)
        for t ∈ T
            i = intervals(f)
            i = i[0.0 .< i .< t]
            @test integral(f, t) ≈ quadgk(f, 0.0, i..., t)[1] rtol = 1e-8
        end

        # B-value
        b = quadgk(t -> integral(f, t)^2, 0.0, Δ + δ)[1]
        @test bvalue_no_q(f) ≈ b rtol = 1e-8
    end

    @testset "DoublePGSE" begin
        δ, Δ, p = 500.0, 2000.0, 1000.0
        tmid = p + δ + Δ
        f = DoublePGSE(δ, Δ, p)

        # Properties
        @test f.δ ≈ δ
        @test f.Δ ≈ Δ
        @test f.tpause ≈ p

        # Function call
        @test f(-100.0) ≈ 0.0
        @test f(prevfloat(0.0)) ≈ 0.0
        @test f(0.0) ≈ 1.0
        @test f(nextfloat(0.0)) ≈ 1.0
        @test f(prevfloat(δ)) ≈ 1.0
        @test f(δ) ≈ 0.0
        @test f(nextfloat(δ)) ≈ 0.0
        @test f(prevfloat(Δ)) ≈ 0.0
        @test f(Δ) ≈ -1.0
        @test f(nextfloat(Δ)) ≈ -1.0
        @test f(prevfloat(Δ + δ)) ≈ -1.0
        @test f(nextfloat(Δ + δ)) ≈ 0.0
        @test f(-prevfloat(tmid)) ≈ 0.0
        @test f(tmid) ≈ 1.0
        @test f(nextfloat(tmid)) ≈ 1.0
        @test f(prevfloat(tmid + δ)) ≈ 1.0
        @test f(nextfloat(tmid + δ)) ≈ 0.0
        @test f(tmid + δ) ≈ 0.0
        @test f(prevfloat(tmid + Δ)) ≈ 0.0
        @test f(tmid + Δ) ≈ -1.0
        @test f(nextfloat(tmid + Δ)) ≈ -1.0
        @test f(prevfloat(tmid + Δ + δ)) ≈ -1.0
        @test f(tmid + Δ + δ) ≈ 0.0
        @test f(nextfloat(tmid + Δ + δ)) ≈ 0.0
        @test f(tmid + Δ + δ + 100.0) ≈ 0.0

        # Integral
        T = LinRange(-δ, p + 2Δ + 3δ, 10)
        for t ∈ T
            i = intervals(f)
            i = i[0.0 .< i .< t]
            @test integral(f, t) ≈ quadgk(f, 0.0, i..., t)[1] rtol = 1e-8
        end

        # B-value
        b = quadgk(t -> integral(f, t)^2, 0.0, tmid + Δ + δ)[1]
        @test bvalue_no_q(f) ≈ b rtol = 1e-8
    end

    @testset "CosOGSE" begin

    end

    @testset "SinOGSE" begin

    end

    @testset "UnknownSequence" begin
        struct UnknownSequence{T} <: TimeProfile{T} end

        f = UnknownSequence{Float64}()

        @test_throws Exception f(0.0)
        @test_throws Exception integral(f, 0.0)
        @test_throws Exception bvalue_no_q(f)
        @test_throws Exception intervals(f)
        @test_throws Exception echotime(f)
        @test_throws Exception get_interval(f, 0.0)
        @test_throws Exception constant_intervals(f)
        @test_throws Exception isconstant(f, 0.0)
    end
end
