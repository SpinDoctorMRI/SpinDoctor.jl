"""
    output_type(problem)

Get type returned by `solve(problem)`.
"""
function output_type end

output_type(::GeneralBTPDE{T}) where {T} = Vector{Complex{T}}
output_type(::IntervalConstantBTPDE{T}) where {T} = Vector{Complex{T}}
output_type(::HADC{T}) where {T} = Vector{T}
output_type(::Karger{T}) where {T} = Complex{T}
output_type(::Laplace{T}) where {T} = NamedTuple
output_type(::MatrixFormalism{T}) where {T} = Vector{Complex{T}}
output_type(::AnalyticalLaplace{T}) where {T} = NamedTuple
output_type(::AnalyticalMatrixFormalism{T}) where {T} = Complex{T}
