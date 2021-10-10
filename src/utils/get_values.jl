"""
    get_values(gradient)

Get q-values and b-values
"""
function get_values(gradient)
    nsequence = length(gradient.sequences)
    if gradient.values_type == "q"
        qvalues = repeat(gradient.values, 1, nsequence)
        bvalues =
            gradient.values .^ 2 .* bvalue_no_q.(gradient.sequences)'
    else
        bvalues = repeat(gradient.values, 1, nsequence)
        qvalues = .√(gradient.values ./ bvalue_no_q.(gradient.sequences)')
    end

    qvalues, bvalues
end
