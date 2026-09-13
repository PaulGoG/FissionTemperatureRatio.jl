# Benchmarks for the parts whose cost scales: the exhaustive breakpoint search, and the level
# density parameter ratio over a full fragmentation range.
#
#     julia --project=bench bench/benchmarks.jl

include(joinpath(@__DIR__, "activate.jl"))

using BenchmarkTools
using FissionTemperatureRatio
using StableRNGs

const MASSES = read_mass_excess(
    joinpath(pkgdir(FissionTemperatureRatio), "data", "mass_excess", "AME2020.ANA")
)
const CHARGES = ChargeDistributionData(
    Dict{Int,Float64}(), Dict{Int,Float64}(), -0.5, 0.6, nothing
)
const PRESCRIPTION = BackShiftedFermiGas(MASSES)

function ratio_sample()
    rng = StableRNG(20260912)
    A_H = collect(126:174)
    value = @. ifelse(A_H ≤ 130, 0.5 - 0.03 * (A_H - 126), 0.38 + 0.0125 * (A_H - 130))
    return (A_H, value .+ 0.004 .* randn(rng, length(A_H)), fill(0.004, length(A_H)))
end

suite = BenchmarkGroup()

let (x, y, σ) = ratio_sample()
    suite["segments"] = BenchmarkGroup()
    for order in 2:5
        suite["segments"]["max_segments=$(order)"] = @benchmarkable fit_segments(
            $x, $y, $σ; max_segments = $order
        )
    end
end

suite["fragmentation"] = @benchmarkable fragmentation_domain(252, 98, 126:174, 5, $CHARGES)

let domain = fragmentation_domain(252, 98, 126:174, 5, CHARGES)
    suite["level_density_ratio"] = BenchmarkGroup()
    for averaging in (RatioOfMeans(), MeanOfRatios())
        suite["level_density_ratio"][string(nameof(typeof(averaging)))] = @benchmarkable(
            level_density_ratio($averaging, $PRESCRIPTION, 252, 98, $domain)
        )
    end
end

tune!(suite)
results = run(suite; verbose = true)
display(median(results))
