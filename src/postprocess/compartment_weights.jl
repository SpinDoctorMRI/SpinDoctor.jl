"""
    compartment_weights(model; TE = 0)

Relative contributions of the compartments to the signal at `b = 0`: the initial spin
density times the volume, attenuated by T₂ relaxation at the echo time `TE`, and normalized
to sum to one.

The ADC of the whole sample in the low-`b` limit is the weighted average of the compartment
ADCs without exchange, e.g. `compartment_weights(model; TE)' * adcs` for compartment ADCs
from [`solve_hadc`](@ref).
"""
function compartment_weights(model; TE = 0)
    volumes = get_cmpt_volumes(model.mesh)
    w = abs.(model.ρ) .* volumes .* exp.(-TE ./ model.T₂)
    w / sum(w)
end
