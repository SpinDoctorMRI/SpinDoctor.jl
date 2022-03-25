"""
    isconstant(profile)

Return `true` if the time profile `profile` is intervalwise constant.
"""
isconstant(::TimeProfile) = false
isconstant(::PGSE) = true
isconstant(::DoublePGSE) = true
