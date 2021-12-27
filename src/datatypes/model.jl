"""
    Model(mesh, ρ, D, T₂, κ, γ)

Discretized biological model with initial spin densities `ρ`, diffusion tensors `D`,
T₂-relaxation times `T₂`, wall permeabilities `κ`, and gyromacnetic ratio `γ` (defaults to
the one of water protons).
"""
Base.@kwdef struct Model{T}
    mesh::FEMesh
    ρ::Vector{Complex{T}}
    D::Vector{SMatrix{3,3,T,9}}
    T₂::Vector{T}
    κ::Vector{T}
    γ::T
end
