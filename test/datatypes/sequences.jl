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
        @test int_F²(f) ≈ b rtol = 1e-8
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
        @test int_F²(f) ≈ b rtol = 1e-8
    end

    @testset "OGSE" begin
        δ, Δ = 2000.0, 3000.0
        for f ∈ (CosOGSE(δ, Δ, 2), SinOGSE(δ, Δ, 2))
            for t ∈ LinRange(0.0, Δ + δ, 7)
                i = filter(s -> 0 < s < t, intervals(f))
                @test integral(f, t) ≈ quadgk(f, 0.0, i..., t; rtol = 1e-12)[1] atol =
                    1e-8 * δ
            end
            b = quadgk(t -> integral(f, t)^2, intervals(f)...; rtol = 1e-12)[1]
            @test int_F²(f) ≈ b rtol = 1e-8
        end
    end

    @testset "Custom profile" begin
        # A profile without closed forms falls back on numerical integration
        struct CustomPGSE <: SpinDoctor.AbstractTimeProfile{Float64}
            pgse::PGSE{Float64}
        end
        (f::CustomPGSE)(t) = f.pgse(t)
        SpinDoctor.intervals(f::CustomPGSE) = intervals(f.pgse)

        pgse = PGSE(500.0, 2000.0)
        f = CustomPGSE(pgse)
        @test echotime(f) ≈ echotime(pgse)
        for t ∈ (300.0, 1000.0, 2300.0)
            @test integral(f, t) isa Float64
            @test integral(f, t) ≈ integral(pgse, t) rtol = 1e-6
        end
        @test int_F²(f) isa Float64
        @test int_F²(f) ≈ int_F²(pgse) rtol = 1e-6
    end
end
