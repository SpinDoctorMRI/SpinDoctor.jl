# Compare two sets of regression outputs written by `record.jl`.
#
# Usage: julia compare.jl <new outputs.toml> <reference outputs.toml>
#
# Prints the relative difference `‖new - ref‖ / ‖ref‖` of every recorded
# quantity and exits with an error if any of them exceeds its tolerance.

using TOML

new = TOML.parsefile(ARGS[1])
ref = TOML.parsefile(ARGS[2])

# Direct computations (assembly, sparse direct solves, closed forms) must agree
# to round-off. Adaptive time stepping and iterative eigensolvers agree up to
# their tolerances.
function tolerance(key)
    if startswith(key, "btpde_qndf") ||
       startswith(key, "btpde_general") ||
       startswith(key, "hadc") ||
       startswith(key, "karger")
        5e-3
    elseif startswith(key, "laplace_mf")
        1e-6
    else
        1e-10
    end
end

flatten(x::Number) = [Float64(x)]
flatten(x::AbstractVector) = Float64.(x)
flatten(x) = error("Cannot flatten $(typeof(x))")

# Relative difference, with a lower bound `scale` on the norm of the reference
function relerr(a, b, scale = 0.0)
    x, y = flatten(a), flatten(b)
    length(x) == length(y) || return Inf
    sqrt(sum(abs2, x - y)) / max(sqrt(sum(abs2, y)), scale, eps())
end

function compare(key, a, b; scale = 0.0)
    if b isa AbstractString
        # The reference failed: report what the new code does instead
        return (a isa AbstractString ? "still fails" : "now runs (reference failed)"), true
    elseif a isa AbstractString
        return "FAILS: " * a, false
    elseif b isa AbstractDict
        msgs = String[]
        ok = true
        ks = sort(collect(keys(b)))
        if "re" in ks && "im" in ks
            # Compare complex signals as a whole: the imaginary parts are often
            # round-off
            haskey(a, "re") && haskey(a, "im") || return "re/im missing", false
            err = relerr([a["re"]; a["im"]], [b["re"]; b["im"]])
            push!(msgs, "signal: $(round(err; sigdigits = 2))")
            ok &= err ≤ tolerance(key)
            filter!(k -> k ∉ ("re", "im"), ks)
        end
        for k in ks
            haskey(a, k) || (push!(msgs, "$k missing"); ok = false; continue)
            msg, o = compare("$key.$k", a[k], b[k])
            push!(msgs, "$k: $msg")
            ok &= o
        end
        return join(msgs, ", "), ok
    else
        err = relerr(a, b, scale)
        isinf(err) && return "different lengths", false
        return string(round(err; sigdigits = 2)), err ≤ tolerance(key)
    end
end

allok = true
for case in sort(collect(keys(ref)))
    ref[case] isa AbstractDict || continue
    println(case)
    if !haskey(new, case)
        println("  missing")
        global allok = false
        continue
    end
    for key in sort(collect(keys(ref[case])))
        if !haskey(new[case], key)
            if ref[case][key] isa AbstractString
                # Error recorded by the reference, not by the new code
                println("  ", rpad(key, 22), "now runs (reference failed)")
            else
                println("  ", rpad(key, 22), "missing")
                global allok = false
            end
            continue
        end
        # Matrix sums can vanish (stiffness matrix, moments of centered
        # geometries): compare them on the scale of the quadratic form
        m = match(r"^(matrix_\w+)_sum$", key)
        scale = isnothing(m) ? 0.0 : abs(get(ref[case], m[1] * "_vAv", 0.0))
        msg, ok = compare(key, new[case][key], ref[case][key]; scale)
        println("  ", rpad(key, 22), ok ? "" : "✗ ", msg)
        global allok &= ok
    end
    for key in sort(collect(setdiff(keys(new[case]), keys(ref[case]))))
        println("  ", rpad(key, 22), "new (not in the reference)")
    end
end
allok || error("Regression outputs differ from the reference")
println("All outputs agree with the reference")
