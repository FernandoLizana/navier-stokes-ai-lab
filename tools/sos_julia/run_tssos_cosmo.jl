# TSSOS + COSMO (chordal decompose) for larger bands where Clarabel OOMs.
# Usage: julia -t auto --project=tssos_env run_tssos_cosmo.jl <dir> [order] [TS] [merge]
using LinearAlgebra
using JSON, NPZ, DynamicPolynomials, TSSOS, JuMP, COSMO

# Use all logical cores for BLAS (PSD eig / projections).
const _NTHR = max(1, something(tryparse(Int, get(ENV, "SOS_NUM_THREADS", "")), Sys.CPU_THREADS))
BLAS.set_num_threads(_NTHR)
# Keep a large working set in RAM (pair with julia --heap-size-hint=…).
# Default OFF: pin_ram caused ReadOnlyMemoryError on D=500 block structure.
const _PIN_RAM = get(ENV, "SOS_PIN_RAM", "0") != "0"
println("threads: Julia=$(Threads.nthreads()) BLAS=$(BLAS.get_num_threads()) CPU=$(Sys.CPU_THREADS) pin_ram=$(_PIN_RAM)"); flush(stdout)

function load_meta(dir); JSON.parsefile(joinpath(dir, "meta.json")); end
function load_A_coo(dir)
    z = npzread(joinpath(dir, "A_coo.npz"))
    D = Int(z["D"]); i = Int.(z["i"]) .+ 1; j = Int.(z["j"]) .+ 1; k = Int.(z["k"]) .+ 1
    return D, i, j, k, Float64.(z["v"])
end
function build_cubic(D, i, j, k, v)
    @polyvar x[1:D]
    mons = [x[i[t]] * x[j[t]] * x[k[t]] for t in eachindex(v)]
    return x, polynomial(v, mons)
end

function main(args)
    isempty(args) && error("usage: run_tssos_cosmo.jl <dir> [order] [TS] [merge]")
    dir = args[1]
    order = length(args) >= 2 ? parse(Int, args[2]) : 2
    TS = length(args) >= 3 ? args[3] : "MD"
    merge = length(args) >= 4 ? (lowercase(args[4]) in ("1","true","yes","merge")) : false
    meta = load_meta(dir)
    println("meta D=$(meta["D"]) C_fullsym=$(get(meta,"C_fullsym",nothing))"); flush(stdout)
    D, i, j, k, v = load_A_coo(dir)
    println("COO nnz=$(length(v)); building poly..."); flush(stdout)
    x, f = build_cubic(D, i, j, k, v)
    h = sum(xi^2 for xi in x) - 1
    pop = [f, h]
    model = Model(optimizer_with_attributes(COSMO.Optimizer,
        "decompose" => true,
        "max_iter" => 10000,
        "eps_abs" => 1e-5,
        "eps_rel" => 1e-5,
        "verbose" => true,
    ))
    println("tssos+COSMO order=$order TS=$TS merge=$merge D=$D ..."); flush(stdout)
    if _PIN_RAM
        # Hold allocations during the heavy SDP; re-enable after.
        GC.enable(false)
    end
    local opt, sol, data
    try
        opt, sol, data = tssos(pop, x, order; numeq=1, TS=TS, merge=merge, QUIET=false, solve=true, model=model)
    finally
        if _PIN_RAM
            GC.enable(true)
            GC.gc(true)
        end
    end
    ub_raw = abs(Float64(opt))
    C_ub = 2 * sqrt(2) * ub_raw
    Cd = Float64(meta["C_dagger"])
    C_fs = get(meta, "C_fullsym", nothing)
    out = Dict(
        "dir" => dir, "D" => D,
        "method" => "TSSOS+COSMO(decompose) order=$order TS=$TS merge=$merge sphere eq",
        "order" => order, "TS" => TS, "merge" => merge,
        "opt" => Float64(opt), "ub_raw_abs" => ub_raw, "C_ub" => C_ub,
        "C_dagger" => Cd, "C_fullsym" => C_fs, "beta" => meta["beta"],
        "closes_vs_Cdagger" => C_ub < Cd - 1e-6,
        "beats_fullsym" => (C_fs !== nothing) && isfinite(C_ub) && (C_ub < Float64(C_fs) - 1e-6),
        "evidence" => "N2_float_SDP", "not_a_theorem" => true, "clay_implication" => "None",
    )
    outfile = joinpath(dir, "cosmo_tssos_result.json")
    open(outfile, "w") do io; JSON.print(io, out, 2); end
    println("FINAL ", JSON.json(out)); println("wrote ", outfile)
end
main(ARGS)
