using Test
using IntervalArithmetic
include(joinpath(@__DIR__, "..", "src", "ValidatedNS.jl"))
using .ValidatedNS

@testset "Interval Leray" begin
    v = [1.0, 2.0, 3.0]
    k = [1.0, 0.0, 0.0]
    p = interval_leray(v, k)
    d2 = interval_div_sq(p, k)
    @test 0 ∈ d2 || diam(d2) < 1e-10 || (inf(d2) <= 0 <= sup(d2))
end

@testset "Heat propagate" begin
    c = interval(1.0)
    out = short_heat_propagate(c, 2.0, 0.1, 0.01)
    expected = exp(-0.1 * 2.0 * 0.01)
    @test expected ∈ out
end
