"""
    solve_multigrad(problem, gradients, args...; kwargs...)

Solve `problem` for multiple magnetic field gradients.
"""
function solve_multigrad(problem, gradients, args...; kwargs...)
    itertimes = zeros(0)
    results = Vector{output_type(problem)}()
    @info "Solving" typeof(problem)
    for grad ∈ gradients
        @info grad
        starttime = time()
        push!(results, solve(problem, grad, args...; kwargs...))
        push!(itertimes, time() - starttime)
    end
    results, itertimes
end
