@with_kw struct Experiment
    gradient::NamedTuple
    btpde::Union{NamedTuple,Nothing} = nothing
    btpde_midpoint::Union{NamedTuple,Nothing} = nothing
    hadc::Union{NamedTuple,Nothing} = nothing
    mf::Union{NamedTuple,Nothing} = nothing
    analytical::Union{NamedTuple,Nothing} = nothing
    karger::Union{NamedTuple,Nothing} = nothing
end
