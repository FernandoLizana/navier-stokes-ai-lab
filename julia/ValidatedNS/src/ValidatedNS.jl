"""ValidatedNS: finite-mode interval probes for NS-MRL Sprint 1.

Evidence: N5 for finite Galerkin identities only — NOT continuum NS.
"""
module ValidatedNS

using IntervalArithmetic
using LinearAlgebra

export interval_leray, interval_div_sq, short_heat_propagate

"""Interval Leray projector for a single real mode k ≠ 0."""
function interval_leray(v::Vector, k::Vector)
    @assert length(v) == 3 && length(k) == 3
    k2 = sum(k .* k)
    if k2 == 0
        return interval.(v)
    end
    vv = interval.(v)
    kk = interval.(k)
    factor = (dot(kk, vv)) / interval(k2)
    return vv .- kk .* factor
end

function interval_div_sq(v_int, k::Vector)
    kk = interval.(k)
    d = dot(kk, v_int)
    return d * d
end

function short_heat_propagate(c::Interval, k2::Real, nu::Real, dt::Real)
    factor = exp(-nu * k2 * dt)
    return c * interval(factor)
end

end # module
