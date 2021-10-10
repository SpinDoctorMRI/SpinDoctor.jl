Base.@kwdef struct Gradient{T}
    directions::Matrix{T}
    sequences::Vector{TimeProfile}
    values::Vector{T}
    values_type::String
end
