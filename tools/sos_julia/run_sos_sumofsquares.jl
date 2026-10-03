# SumOfSquares upper bound on max A(z,z,z) over the sphere via ideal membership:
#   γ - f - λ(x)*(||x||^2 - 1) ∈ SOSCone
# with λ ∈ ℝ[x]_1. (Odd cubic cancelled by λ·quadratic → even SOS remainder.)
#
# Usage: julia --project=. run_sos_sumofsquares.jl data/band_123
# FINITE Galerkin only. Float SDP = N2. Not Clay.

using JSON
using NPZ
using DynamicPolynomials
using SumOfSquares
using Clarabel

function load_meta(dir::String)
    JSON.parsefile(joinpath(dir, "meta.json"))
end

function load_A_coo(dir::String)
    z = npzread(joinpath(dir, "A_coo.npz"))
    D = Int(z["D"])
    i = Int.(z["i"]) .+ 1
    j = Int.(z["j"]) .+ 1
    k = Int.(z["k"]) .+ 1
    v = Float64.(z["v"])
    return D, i, j, k, v
end

function build_cubic(D, i, j, k, v)
    @polyvar x[1:D]
    f = zero(x[1] * x[1] * x[1])
    @inbounds for t in eachindex(v)
        f += v[t] * x[i[t]] * x[j[t]] * x[k[t]]
    end
    return x, f
end

function sos_ub(x, f; quiet::Bool=true)
    model = SOSModel(Clarabel.Optimizer)
    if quiet
        set_silent(model)
    end
    @variable(model, γ)
    λ = @variable(model, variable_type = Poly(monomials(x, 0:1)))
    S = sum(xi^2 for xi in x)
    # Ideal reduction: γ - f ≡ SOS  (mod ||x||^2-1)
    @constraint(model, γ - f - λ * (S - 1) in SOSCone())
    @objective(model, Min, γ)
    optimize!(model)
    return objective_value(model), termination_status(model), primal_status(model), dual_status(model)
end

function main(args)
    isempty(args) && error("usage: julia --project=. run_sos_sumofsquares.jl <data_dir>")
    dir = args[1]
    meta = load_meta(dir)
    println("meta D=$(meta["D"]) C_fullsym=$(get(meta,"C_fullsym",nothing)) beta=$(meta["beta"])")
    flush(stdout)
    D, i, j, k, v = load_A_coo(dir)
    println("COO nnz=$(length(v)); building poly...")
    flush(stdout)
    x, f = build_cubic(D, i, j, k, v)
    println("SOS ub for +f ..."); flush(stdout)
    ub_pos, st_pos, pr_pos, du_pos = sos_ub(x, f)
    println("  ub_pos=$(ub_pos) term=$(st_pos) prim=$(pr_pos) dual=$(du_pos)"); flush(stdout)
    println("SOS ub for -f ..."); flush(stdout)
    ub_neg, st_neg, pr_neg, du_neg = sos_ub(x, -f)
    println("  ub_neg=$(ub_neg) term=$(st_neg) prim=$(pr_neg) dual=$(du_neg)"); flush(stdout)
    ub_raw = max(abs(Float64(ub_pos)), abs(Float64(ub_neg)))
    Cd = Float64(meta["C_dagger"])
    C_ub = 2 * sqrt(2) * ub_raw
    C_fs = get(meta, "C_fullsym", nothing)
    out = Dict(
        "dir" => dir,
        "D" => D,
        "method" => "SumOfSquares+Clarabel γ-f-λ(||x||^2-1) in SOSCone",
        "ub_pos" => Float64(ub_pos),
        "ub_neg" => Float64(ub_neg),
        "status_pos" => string(st_pos),
        "status_neg" => string(st_neg),
        "ub_raw_abs" => ub_raw,
        "C_ub" => C_ub,
        "C_dagger" => Cd,
        "C_fullsym" => C_fs,
        "beta" => meta["beta"],
        "closes_vs_Cdagger" => C_ub < Cd - 1e-6,
        "beats_fullsym" => (C_fs !== nothing) && (C_ub < Float64(C_fs) - 1e-6),
        "evidence" => "N2_float_SDP",
        "not_a_theorem" => true,
        "clay_implication" => "None",
    )
    outfile = joinpath(dir, "sos_sos_result.json")
    open(outfile, "w") do io
        JSON.print(io, out, 2)
    end
    println("FINAL ", JSON.json(out))
    println("wrote ", outfile)
end

main(ARGS)
