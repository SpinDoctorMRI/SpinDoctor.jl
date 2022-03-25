"""
    compute_signal(M, ξ)

Compute signal from magnetization `ξ`, using the mass matrix `M` for integration.

Given a mesh `mesh` and a vector of compartment mass matrices `M_cmpts`, a vector of
compartment signals may be obtained by `compute_signal.(M_cmpts, split_field(mesh, ξ))`.
"""
compute_signal(M, ξ) = sum(M * ξ)
