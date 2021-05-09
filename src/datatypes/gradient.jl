@with_kw struct Gradient
    directions::Array
    sequences::Vector{TimeProfile}
    values::Vector{Float64}
    values_type::String
end