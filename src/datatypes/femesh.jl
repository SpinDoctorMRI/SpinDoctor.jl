struct FEMesh{T,dim}
    point_map::Vector{Vector{Int}}
    points::Vector{Matrix{T}}
    facets::Matrix{Matrix{Int}}
    elements::Vector{Matrix{Int}}
end

function FEMesh(; point_map, points, facets, elements)
    T = eltype(first(points))
    dim = size(first(points), 1)
    FEMesh{T,dim}(point_map, points, facets, elements)
end
