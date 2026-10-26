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
        1e-3
    elseif startswith(key, "laplace_mf")
        1e-6
    else
        1e-10
    end
end

flatten(x::Number) = [Float64(x)]
flatten(x::AbstractVector) = Float64.(x)
flatten(x) = error("Cannot flatten $(typeof(x))")

function compare(key, a, b)
    if b isa AbstractString
        # The reference failed: report what the new code does instead
        return (a isa AbstractString ? "still fails" : "now runs (reference failed)"), true
    elseif a isa AbstractString
        return "FAILS: " * a, false
    elseif b isa AbstractDict
        msgs = String[]
        ok = true
        for k in sort(collect(keys(b)))
            haskey(a, k) || (push!(msgs, "$k missing"); ok = false; continue)
            msg, o = compare("$key.$k", a[k], b[k])
            push!(msgs, "$k: $msg")
            ok &= o
        end
        return join(msgs, ", "), ok
    else
        x, y = flatten(a), flatten(b)
        length(x) == length(y) || return "length $(length(x)) ≠ $(length(y))", false
        err = sqrt(sum(abs2, x - y)) / max(sqrt(sum(abs2, y)), eps())
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
            println("  ", rpad(key, 22), "missing")
            global allok = false
            continue
        end
        msg, ok = compare(key, new[case][key], ref[case][key])
        println("  ", rpad(key, 22), ok ? "" : "✗ ", msg)
        global allok &= ok
    end
end
allok || error("Regression outputs differ from the reference")
println("All outputs agree with the reference")
