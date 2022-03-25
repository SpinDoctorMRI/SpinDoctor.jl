"""
    initial_conditions(model)

Get initial conditions at all degrees of freedom.
"""
initial_conditions(model) =
    mapreduce((ρ, p) -> fill(ρ, size(p, 2)), vcat, model.ρ, model.mesh.points)
