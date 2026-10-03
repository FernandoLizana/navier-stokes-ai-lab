using Pkg

# Core SDP / SOS stack from General
Pkg.add(["Clarabel", "DynamicPolynomials", "JSON", "JuMP", "NPZ", "SumOfSquares", "MutableArithmetics"])

# TSSOS is not always in General — install from upstream
try
    Pkg.add(url="https://github.com/wangjie212/TSSOS")
catch e
    @warn "TSSOS GitHub add failed" exception=e
end

try
    Pkg.add(url="https://github.com/wangjie212/NCTSSOS")
catch e
    @warn "NCTSSOS GitHub add failed" exception=e
end

Pkg.status()
