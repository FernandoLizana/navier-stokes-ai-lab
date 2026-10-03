using NPZ, DynamicPolynomials
z = npzread("data/band_1to12/A_coo.npz")
D = Int(z["D"]); i = Int.(z["i"]) .+ 1; j = Int.(z["j"]) .+ 1; k = Int.(z["k"]) .+ 1; v = Float64.(z["v"])
println("D=$D nnz=$(length(v))")
@polyvar x[1:D]
t0 = time()
mons = [x[i[t]] * x[j[t]] * x[k[t]] for t in eachindex(v)]
println("mons $(round(time()-t0;digits=2))s")
t0 = time()
f = polynomial(v, mons)
println("poly $(round(time()-t0;digits=2))s nterms=$(nterms(f))")
