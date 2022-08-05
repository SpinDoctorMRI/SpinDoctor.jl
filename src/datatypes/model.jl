"""
    Model(; mesh, ρ, D, T₂, κ, γ)

Finite element discretized biological model with initial spin densities `ρ`,
diffusion tensors `D`, T₂-relaxation times `T₂`, wall permeabilities `κ`, and
gyromagnetic ratio `γ`. The vectors `ρ`, `D`, `T₂`, are of length
`ncompartment`, while `κ` is of length `nboundary`.
"""
Base.@kwdef struct Model{T,dim}
    mesh::FEMesh{T,dim}
    ρ::Vector{Complex{T}}
    D::Vector{Matrix{T}}
    T₂::Vector{T}
    κ::Vector{T}
    γ::T
end
