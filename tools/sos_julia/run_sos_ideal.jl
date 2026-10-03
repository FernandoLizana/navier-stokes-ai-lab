# Ideal-based SOS: γ - f ∈ SOSCone on the quotient ℝ[x]/(||x||^2-1)
using JSON
using NPZ
using DynamicPolynomials
using SumOfSquares
using Clarabel
using SemialgebraicSets

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
    S = sum(xi^2 for xi in x)
    I = Ideal(S - 1)
    @constraint(model, γ - f in SOSCone(), domain = I)
    @objective(model, Min, γ)
    optimize!(model)
    return objective_value(model), termination_status(model), primal_status(model)
end

function main(args)
    dir = args[1]
    meta = load_meta(dir)
    println("meta D=$(meta["D"]) C_fullsym=$(meta["C_fullsym"]) beta=$(meta["beta"])"); flush(stdout)
    D, i, j, k, v = load_A_coo(dir)
    x, f = build_cubic(D, i, j, k, v)
    println("SOS Ideal +f ..."); flush(stdout)
    ub_pos, st_pos, pr_pos = sos_ub(x, f)
    println("  ub_pos=$(ub_pos) $(st_pos) $(pr_pos)"); flush(stdout)
    println("SOS Ideal -f ..."); flush(stdout)
    ub_neg, st_neg, pr_neg = sos_ub(x, -f)
    println("  ub_neg=$(ub_neg) $(st_neg) $(pr_neg)"); flush(stdout)
    ub_raw = max(abs(Float64(ub_pos)), abs(Float64(ub_neg)))
    Cd = Float64(meta["C_dagger"])
    C_ub = 2 * sqrt(2) * ub_raw
    out = Dict(
        "D" => D,
        "method" => "SumOfSquares Ideal(||x||^2-1) SOSCone",
        "ub_pos" => Float64(ub_pos),
        "ub_neg" => Float64(ub_neg),
        "status_pos" => string(st_pos),
        "status_neg" => string(st_neg),
        "ub_raw_abs" => ub_raw,
        "C_ub" => C_ub,
        "C_dagger" => Cd,
        "C_fullsym" => meta["C_fullsym"],
        "beta" => meta["beta"],
        "closes_vs_Cdagger" => C_ub < Cd - 1e-6,
        "evidence" => "N2_float_SDP",
        "not_a_theorem" => true,
    )
    outfile = joinpath(dir, "sos_ideal_result.json")
    open(outfile, "w") do io
        JSON.print(io, out, 2)
    end
    println("FINAL ", JSON.json(out))
end

main(ARGS)
