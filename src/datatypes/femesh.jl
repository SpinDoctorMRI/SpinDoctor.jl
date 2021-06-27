@with_kw struct FEMesh{T}
    point_map::Vector{Vector{Int}}
    points::Vector{Matrix{T}}
    facets::Matrix{Matrix{Int}}
    elements::Vector{Matrix{Int}}
end
