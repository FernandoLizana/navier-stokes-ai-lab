# tools/sos_julia/run_sos_band.jl
# TSSOS on physical cubic A(z,z,z) with ‖z‖²=1.
# Need ub < beta = C_†/(2√2) ≈ 3.3807 for all-IC low-slab via L-0027.
#
# Usage:
#   julia --project=. run_sos_band.jl data/band_123456
#
# FINITE Galerkin only. Float SDP = numerical evidence, not a theorem. Not Clay.

using JSON
using NPZ
using DynamicPolynomials
using TSSOS
using JuMP
using Clarabel

function load_meta(dir::String)
    JSON.parsefile(joinpath(dir, "meta.json"))
end

function load_A_coo(dir::String)
    z = npzread(joinpath(dir, "A_coo.npz"))
    D = Int(z["D"])
    i = Int.(z["i"]) .+ 1  # Julia 1-based
    j = Int.(z["j"]) .+ 1
    k = Int.(z["k"]) .+ 1
    v = Float64.(z["v"])
    return D, i, j, k, v
end

function build_cubic(D, i, j, k, v)
    @polyvar x[1:D]
    f = zero(x[1])
    for t in eachindex(v)
        f += v[t] * x[i[t]] * x[j[t]] * x[k[t]]
    end
    return x, f
end

function run_tssos(x, f; order::Int=2, quiet::Bool=false)
    # Pop: max f s.t. ||x||^2 = 1  ↔  equality g=1-sum(x.^2)
    g = [1 - sum(x .^ 2)]
    # Sign-symmetric: solve max f and max -f
    results = Dict{String,Any}()
    for (name, poly) in (("pos", f), ("neg", -f))
        println("TSSOS $(name) order=$(order) D=$(length(x)) ...")
        flush(stdout)
        opt, sol, data = tssos_first(
            poly,
            x,
            order;
            numeq=1,
            equality=g,
            TS="block",
            QUIET=quiet,
            solve=true,
            Motzkin=false,
        )
        results[name] = Dict("opt" => opt, "sol_status" => string(something(sol, "none")))
        println("  $(name) opt=$(opt)")
        flush(stdout)
    end
    return results
end

function main(args)
    isempty(args) && error("usage: julia --project=. run_sos_band.jl <data_dir>")
    dir = args[1]
    meta = load_meta(dir)
    println("meta: ", JSON.json(meta))
    flush(stdout)
    D, i, j, k, v = load_A_coo(dir)
    println("loaded A COO D=$(D) nnz=$(length(v))")
    flush(stdout)
    x, f = build_cubic(D, i, j, k, v)
    order = length(args) >= 2 ? parse(Int, args[2]) : 2
    res = run_tssos(x, f; order=order, quiet=false)
    beta = Float64(meta["beta"])
    Cd = Float64(meta["C_dagger"])
    # best upper on |A| — TSSOS returns lower bounds on the max for the dual?
    # TSSOS opt is typically an upper bound on the POP maximum (relaxation).
    ub_raw = max(abs(Float64(res["pos"]["opt"])), abs(Float64(res["neg"]["opt"])))
    # If solver returns -Inf/NaN handle
    C_ub = 2 * sqrt(2) * ub_raw
    closes = C_ub < Cd - 1e-6
    out = Dict(
        "dir" => dir,
        "D" => D,
        "nnz" => length(v),
        "order" => order,
        "beta" => beta,
        "C_dagger" => Cd,
        "C_fullsym" => get(meta, "C_fullsym", nothing),
        "tssos" => res,
        "ub_raw_abs" => ub_raw,
        "C_ub" => C_ub,
        "closes_band_vs_Cdagger" => closes,
        "evidence" => "N2_float_SDP",
        "not_a_theorem" => true,
        "clay_implication" => "None",
    )
    outfile = joinpath(dir, "sos_result.json")
    open(outfile, "w") do io
        JSON.print(io, out, 2)
    end
    println("FINAL ", JSON.json(out))
    println("wrote ", outfile)
end

main(ARGS)
