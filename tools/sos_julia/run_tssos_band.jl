# TSSOS hierarchy on physical cubic with sphere equality.
# Usage: julia --project=. ../run_tssos_band.jl ../data/band_123
# (run from tssos_env/)
# FINITE Galerkin only. Float SDP = N2. Not Clay.

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
    i = Int.(z["i"]) .+ 1
    j = Int.(z["j"]) .+ 1
    k = Int.(z["k"]) .+ 1
    v = Float64.(z["v"])
    return D, i, j, k, v
end

function build_cubic(D, i, j, k, v)
    @polyvar x[1:D]
    # `polynomial(coeffs, mons)` — Dict/sum paths are too slow for nnz≳1e5.
    mons = [x[i[t]] * x[j[t]] * x[k[t]] for t in eachindex(v)]
    return x, polynomial(v, mons)
end

function run_one(x, f, order; TS="MD", merge=false, md=3, loose=false, dualize=false)
    h = sum(xi^2 for xi in x) - 1
    pop = [f, h]
    println("tssos order=$(order) TS=$(TS) merge=$(merge) md=$(md) dualize=$(dualize) D=$(length(x)) Clarabel ..."); flush(stdout)
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    if loose
        # Coarser but cheaper: enough to compare vs C_fullsym / C_†.
        set_optimizer_attribute(model, "tol_gap_abs", 1e-5)
        set_optimizer_attribute(model, "tol_gap_rel", 1e-5)
        set_optimizer_attribute(model, "tol_feas", 1e-5)
        set_optimizer_attribute(model, "max_iter", 80)
    end
    GC.gc()
    # TSSOS ≥1.5 uses `tssos`; pass Clarabel via model= (default is Mosek+license).
    opt, sol, data = tssos(pop, x, order; numeq=1, TS=TS, merge=merge, md=md, dualize=dualize, QUIET=false, solve=true, model=model)
    println("  opt=$(opt)"); flush(stdout)
    return opt, sol, data
end

function main(args)
    isempty(args) && error("usage: julia --project=tssos_env run_tssos_band.jl <data_dir> [order] [TS] [merge] [loose] [dualize]")
    dir = args[1]
    order = length(args) >= 2 ? parse(Int, args[2]) : 2
    TS = length(args) >= 3 ? args[3] : "MD"
    merge = length(args) >= 4 ? (lowercase(args[4]) in ("1", "true", "yes", "merge")) : false
    loose = length(args) >= 5 ? (lowercase(args[5]) in ("1", "true", "yes", "loose")) : false
    dualize = length(args) >= 6 ? (lowercase(args[6]) in ("1", "true", "yes", "dualize")) : false
    meta = load_meta(dir)
    println("meta D=$(meta["D"]) C_fullsym=$(get(meta,"C_fullsym",nothing)) beta=$(meta["beta"])")
    flush(stdout)
    D, i, j, k, v = load_A_coo(dir)
    println("COO nnz=$(length(v)); building poly..."); flush(stdout)
    x, f = build_cubic(D, i, j, k, v)
    # Odd cubic on centrally symmetric sphere ⇒ max f = −min f; one SDP suffices.
    opt_pos, _, _ = run_one(x, f, order; TS=TS, merge=merge, loose=loose, dualize=dualize)
    opt_neg = opt_pos
    ub_raw = abs(Float64(opt_pos))
    Cd = Float64(meta["C_dagger"])
    C_ub = 2 * sqrt(2) * ub_raw
    C_fs = get(meta, "C_fullsym", nothing)
    out = Dict(
        "dir" => dir,
        "D" => D,
        "method" => "TSSOS tssos order=$(order) TS=$(TS) merge=$(merge) dualize=$(dualize) Clarabel sphere equality",
        "order" => order,
        "TS" => TS,
        "merge" => merge,
        "loose" => loose,
        "dualize" => dualize,
        "opt_pos" => Float64(opt_pos),
        "opt_neg" => Float64(opt_neg),
        "ub_raw_abs" => ub_raw,
        "C_ub" => C_ub,
        "C_dagger" => Cd,
        "C_fullsym" => C_fs,
        "beta" => meta["beta"],
        "closes_vs_Cdagger" => C_ub < Cd - 1e-6,
        "beats_fullsym" => (C_fs !== nothing) && isfinite(C_ub) && (C_ub < Float64(C_fs) - 1e-6),
        "evidence" => "N2_float_SDP",
        "not_a_theorem" => true,
        "clay_implication" => "None",
    )
    outfile = joinpath(dir, "tssos_result.json")
    open(outfile, "w") do io
        JSON.print(io, out, 2)
    end
    println("FINAL ", JSON.json(out))
    println("wrote ", outfile)
end

main(ARGS)
