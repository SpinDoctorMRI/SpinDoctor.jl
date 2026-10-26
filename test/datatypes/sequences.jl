@testset "sequences.jl" begin
    @testset "PGSE" begin
        delta, Delta = 500.0, 2000.0
        f = PGSE(delta, Delta)
        @test f.delta ≈ delta
        @test f.Delta ≈ Delta
        @test f(-100.0) ≈ 0.0
        @test f(prevfloat(0.0)) ≈ 0.0
        @test f(0.0) ≈ 1.0
        @test f(prevfloat(delta)) ≈ 1.0
        @test f(delta) ≈ 0.0
        @test f(nextfloat(delta)) ≈ 0.0
        @test f(prevfloat(Delta)) ≈ 0.0
        @test f(Delta) ≈ -1.0
        @test f(nextfloat(Delta)) ≈ -1.0
        @test f(Delta + delta + 100.0) ≈ 0.0

        # Integral
        T = LinRange(-delta, Delta + 2delta, 10)
        for t in T
            i = intervals(f)
            i = i[0.0 .< i .< t]
            @test integral(f, t) ≈ quadgk(f, 0.0, i..., t)[1] rtol = 1e-8
        end

        # B-value
        b = quadgk(t -> integral(f, t)^2, 0.0, Delta + delta)[1]
        @test int_F2(f) ≈ b rtol = 1e-8
    end

    @testset "DoublePGSE" begin
        delta, Delta, p = 500.0, 2000.0, 1000.0
        tmid = p + delta + Delta
        f = DoublePGSE(delta, Delta, p)

        # Properties
        @test f.delta ≈ delta
        @test f.Delta ≈ Delta
        @test f.tpause ≈ p

        # Function call
        @test f(-100.0) ≈ 0.0
        @test f(prevfloat(0.0)) ≈ 0.0
        @test f(0.0) ≈ 1.0
        @test f(nextfloat(0.0)) ≈ 1.0
        @test f(prevfloat(delta)) ≈ 1.0
        @test f(delta) ≈ 0.0
        @test f(nextfloat(delta)) ≈ 0.0
        @test f(prevfloat(Delta)) ≈ 0.0
        @test f(Delta) ≈ -1.0
        @test f(nextfloat(Delta)) ≈ -1.0
        @test f(prevfloat(Delta + delta)) ≈ -1.0
        @test f(nextfloat(Delta + delta)) ≈ 0.0
        @test f(-prevfloat(tmid)) ≈ 0.0
        @test f(tmid) ≈ 1.0
        @test f(nextfloat(tmid)) ≈ 1.0
        @test f(prevfloat(tmid + delta)) ≈ 1.0
        @test f(nextfloat(tmid + delta)) ≈ 0.0
        @test f(tmid + delta) ≈ 0.0
        @test f(prevfloat(tmid + Delta)) ≈ 0.0
        @test f(tmid + Delta) ≈ -1.0
        @test f(nextfloat(tmid + Delta)) ≈ -1.0
        @test f(prevfloat(tmid + Delta + delta)) ≈ -1.0
        @test f(tmid + Delta + delta) ≈ 0.0
        @test f(nextfloat(tmid + Delta + delta)) ≈ 0.0
        @test f(tmid + Delta + delta + 100.0) ≈ 0.0

        # Integral
        T = LinRange(-delta, p + 2Delta + 3delta, 10)
        for t in T
            i = intervals(f)
            i = i[0.0 .< i .< t]
            @test integral(f, t) ≈ quadgk(f, 0.0, i..., t)[1] rtol = 1e-8
        end

        # B-value
        b = quadgk(t -> integral(f, t)^2, 0.0, tmid + Delta + delta)[1]
        @test int_F2(f) ≈ b rtol = 1e-8
    end

    @testset "OGSE" begin
        delta, Delta = 2000.0, 3000.0
        for f in (CosOGSE(delta, Delta, 2), SinOGSE(delta, Delta, 2))
            for t in LinRange(0.0, Delta + delta, 7)
                i = filter(s -> 0 < s < t, intervals(f))
                @test integral(f, t) ≈ quadgk(f, 0.0, i..., t; rtol = 1e-12)[1] atol =
                    1e-8 * delta
            end
            b = quadgk(t -> integral(f, t)^2, intervals(f)...; rtol = 1e-12)[1]
            @test int_F2(f) ≈ b rtol = 1e-8
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
        for t in (300.0, 1000.0, 2300.0)
            @test integral(f, t) isa Float64
            @test integral(f, t) ≈ integral(pgse, t) rtol = 1e-6
        end
        @test int_F2(f) isa Float64
        @test int_F2(f) ≈ int_F2(pgse) rtol = 1e-6
    end
end
