# Benchmarks for the parts whose cost scales: the exhaustive breakpoint search, and the inversion
# of the relation between R_T and E*_H/TXE over a full fragmentation range.
#
#     julia bench/benchmarks.jl

include(joinpath(@__DIR__, "activate.jl"))

using BenchmarkTools
using FissionFragmentsDomain
using FissionTemperatureRatio
using StableRNGs

const MASSES = read_mass_excess_table(String(AME2020_MASS_EXCESS_FILE))
const SYSTEM = spontaneous_fission(Nuclide(98, 252))
const CHARGES = charge_model(MASSES, SYSTEM)
const MODEL = BackShiftedFermiGas(MASSES)

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

suite["fragmentation"] = @benchmarkable fragmentation_domain($SYSTEM, $CHARGES, 126:174)

# The inversion of a full r_ν(A_H) curve, as a run applies it to every tabulated curve.
let domain = fragmentation_domain(SYSTEM, CHARGES, 126:174)
    A_H = collect(126:174)
    r_ν = RatioCurve(
        A_H,
        @.(ifelse(A_H ≤ 130, 0.5 - 0.03 * (A_H - 126), 0.38 + 0.0125 * (A_H - 130))),
        fill(0.004, length(A_H)),
        "reference",
    )
    energies = Dict(A => 178.0 - 0.02 * (A - 132)^2 for A in A_H)
    weighted = ChargeResolved(mean_total_excitation(MASSES, domain, energies))
    suite["temperature_ratio"] = BenchmarkGroup()
    for (name, averaging) in (
        "charge_resolved" => ChargeResolved(),
        "charge_resolved_weighted" => weighted,
        "ratio_of_means" => RatioOfMeans(),
    )
        suite["temperature_ratio"][name] = @benchmarkable(
            temperature_ratio($averaging, $MODEL, $domain, $r_ν)
        )
    end
end

tune!(suite)
results = run(suite; verbose = true)
display(median(results))
