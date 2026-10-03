# Export TSSOS Gram matrices + rationalize for N5 verification (band_123 D=52).
# FINITE Galerkin only. Not Clay.
using JSON
using NPZ
using DynamicPolynomials
using TSSOS
using JuMP
using Clarabel
using LinearAlgebra

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
    mons = [x[i[t]] * x[j[t]] * x[k[t]] for t in eachindex(v)]
    return x, polynomial(v, mons)
end

function exp_vec_to_list(exps::Vector{UInt16}, D::Int)
    cnt = zeros(Int, D)
    for idx in exps
        cnt[Int(idx)] += 1
    end
    return cnt
end

function rationalize_matrix(Q::AbstractMatrix{<:Real}, denom::Int)
    R = Matrix{BigInt}(undef, size(Q)...)
    err = 0.0
    invd = 1.0 / denom
    for r in axes(Q, 1), c in axes(Q, 2)
        val = Float64(Q[r, c])
        q = round(BigInt, val * denom)
        R[r, c] = q
        err = max(err, abs(val - Float64(q) * invd))
    end
    return R, err
end

function psd_min_eig(Q::AbstractMatrix{<:Real})
    ev = eigvals(Symmetric(Q))
    return minimum(ev)
end

function _jsonify_multiplier(mul)
    if mul isa Real
        return Float64(mul)
    elseif mul isa AbstractVector
        return [_jsonify_multiplier(x) for x in mul]
    else
        return Float64(mul)
    end
end

function _read_tssos_prefs(dir::String)
    merge = false
    loose = false
    path = joinpath(dir, "tssos_result.json")
    if isfile(path)
        d = JSON.parsefile(path)
        merge = get(d, "merge", false)
        loose = get(d, "loose", false)
    end
    return merge, loose
end

function export_gram(
    dir::String;
    order::Int = 2,
    TS::String = "MD",
    denom::Int = 10^8,
    psd_tol::Float64 = 1e-7,
    merge::Union{Nothing,Bool} = nothing,
    loose::Union{Nothing,Bool} = nothing,
)
    meta = load_meta(dir)
    pref_merge, pref_loose = _read_tssos_prefs(dir)
    use_merge = merge === nothing ? pref_merge : merge
    use_loose = loose === nothing ? pref_loose : loose
    D, i, j, k, v = load_A_coo(dir)
    x, f = build_cubic(D, i, j, k, v)
    h = sum(xi^2 for xi in x) - 1
    pop = [f, h]
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    if use_loose
        set_optimizer_attribute(model, "tol_gap_abs", 1e-5)
        set_optimizer_attribute(model, "tol_gap_rel", 1e-5)
        set_optimizer_attribute(model, "tol_feas", 1e-5)
        set_optimizer_attribute(model, "max_iter", 120)
    end
    println("export_gram D=$D merge=$use_merge loose=$use_loose Gram=true ..."); flush(stdout)
    opt, sol, data = tssos(
        pop,
        x,
        order;
        numeq = 1,
        TS = TS,
        merge = use_merge,
        md = 3,
        dualize = false,
        GroebnerBasis = false,
        QUIET = false,
        solve = true,
        Gram = true,
        solution = false,
        model = model,
    )
    ub_raw = abs(Float64(opt))
    Cd = Float64(meta["C_dagger"])
    C_ub = 2 * sqrt(2) * ub_raw

    basis_monomials = Vector{Vector{Int}}(undef, length(data.basis[1]))
    for (bi, exps) in enumerate(data.basis[1])
        basis_monomials[bi] = exp_vec_to_list(exps, D)
    end

    gram_out = Vector{Any}()
    max_rat_err = 0.0
    min_psd = Inf
    if data.GramMat === nothing
        error("GramMat missing; rerun with Gram=true")
    end
    for (ci, gcons) in enumerate(data.GramMat)
        for (ki, Qraw) in enumerate(gcons)
            Q = Qraw isa Matrix ? Qraw : reshape([Float64(Qraw)], 1, 1)
            lam_min = psd_min_eig(Q)
            min_psd = min(min_psd, lam_min)
            R, err = rationalize_matrix(Q, denom)
            max_rat_err = max(max_rat_err, err)
            Rf = Matrix{Float64}(R) ./ denom
            lam_min_r = psd_min_eig(Rf)
            block_idx = ci <= length(data.blocks) ? data.blocks[ci][ki] : Int[]
            push!(
                gram_out,
                Dict(
                    "constraint" => ci,
                    "block_k" => ki,
                    "block_monomial_indices" => block_idx,
                    "size" => size(Q, 1),
                    "denom" => denom,
                    "numerators" => [Int.(R[r, :]) for r in axes(R, 1)],
                    "psd_min_eig_float" => lam_min,
                    "psd_min_eig_rationalized" => lam_min_r,
                ),
            )
        end
    end
    if !isfinite(min_psd)
        error("no Gram blocks exported")
    end

    multiplier = data.multiplier === nothing ? Any[] : _jsonify_multiplier(data.multiplier)
    out = Dict(
        "dir" => dir,
        "D" => D,
        "order" => order,
        "TS" => TS,
        "merge" => use_merge,
        "loose" => use_loose,
        "opt" => Float64(opt),
        "ub_raw_abs" => ub_raw,
        "C_ub" => C_ub,
        "C_dagger" => Cd,
        "C_fullsym" => get(meta, "C_fullsym", nothing),
        "beta" => meta["beta"],
        "denom" => denom,
        "max_rationalization_error" => max_rat_err,
        "min_psd_eig_all_blocks" => min_psd,
        "multiplier" => multiplier,
        "basis_monomials" => basis_monomials,
        "gram_blocks" => gram_out,
        "n_gram_blocks" => length(gram_out),
        "evidence" => "N5_gram_rationalized",
        "clay_implication" => "None",
    )
    outfile = joinpath(dir, "gram_export.json")
    open(outfile, "w") do io
        JSON.print(io, out, 2)
    end
    println("wrote ", outfile)
    println("opt=", opt, " C_ub=", C_ub, " min_psd=", min_psd, " max_rat_err=", max_rat_err)
    if min_psd < -psd_tol
        error("PSD check failed: min eig = $(min_psd)")
    end
    if C_ub >= Cd - 1e-6
        error("C_ub not below C_dagger")
    end
    return out
end

dir = length(ARGS) >= 1 ? ARGS[1] : "data/band_123"
merge_arg = length(ARGS) >= 2 ? (lowercase(ARGS[2]) in ("1", "true", "yes", "merge")) : nothing
loose_arg = length(ARGS) >= 3 ? (lowercase(ARGS[3]) in ("1", "true", "yes", "loose")) : nothing
export_gram(dir; merge=merge_arg, loose=loose_arg)
