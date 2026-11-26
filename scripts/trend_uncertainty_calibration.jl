# Calibration of the uncertainty of the systematic trend in simulated pools:
#
#     julia scripts/trend_uncertainty_calibration.jl <workers> <ar1 | nugget> <replicates> <pool> …
#
# A pool is `<k>w`, k datasets measured at every mass number, or `<k>i`, k datasets measured on a
# grid of 1.5 u and interpolated onto the mass numbers; for example `4w 14w 6i`. The replicates
# are shared among `<workers>` processes, or run in this one for 0; a replicate of four datasets
# takes about half a minute on one core, and the time grows with the number of datasets. The
# method page quotes the tables this prints for `ar1 80 4w 6w 6i 14w` and `nugget 80 4w 6w 6i 14w`.
#
# Every dataset follows one curve of three segments over A_H = 126 … 174, with an error that is
# correlated along the mass axis and a quoted statistical uncertainty of 0.004. The error is
# AR(1) with ρ = 0.8 and a standard deviation of 0.03 (`ar1`), or of 0.026 beside an
# uncorrelated part of 0.015 (`nugget`). A replicate pools k such datasets, fits ρ to the
# correlogram of their deviations, selects the number of segments and fits the trend as a run
# does, and refits it with each dataset left out. Three uncertainties of the fitted r_ν follow:
# the propagated one, the delete-one jackknife over the refits, and the larger of the two, which
# is what a run reports.
#
# Two things are scored over the replicates. For a yield-weighted mean of the curve, which stands
# for ⟨R_T⟩: the standard deviation of the mean over the replicates divided by the root mean
# square of the uncertainty reported for it, 1 for a calibrated uncertainty and above 1 for one
# that is too small. For the curve at each mass number: the root mean square of
# z = (r̂ − r)/σ and the fractions of |z| below 1 and below 2, which are 0.683 and 0.954 for a
# calibrated normal error. The seeds are fixed, so the same arguments print the same table
# whatever the number of workers.
#
# The pooling and the jackknife are the package's own internal functions, so that what is scored
# is what a run computes. Nothing is written.

include(joinpath(@__DIR__, "activate.jl"))

using Distributed: @everywhere, addprocs, pmap
using Printf: @printf, @sprintf
using Statistics: mean, std

const USAGE = "usage: julia scripts/trend_uncertainty_calibration.jl <workers> \
               <ar1 | nugget> <replicates> <pool> …"

# The workers run in this environment, which the activation above has resolved already.
if length(ARGS) ≥ 4 && tryparse(Int, ARGS[1]) !== nothing && parse(Int, ARGS[1]) > 0
    addprocs(parse(Int, ARGS[1]); exeflags = `--project=$(@__DIR__)`)
end

@everywhere using FissionTemperatureRatio:
    FissionTemperatureRatio,
    RatioCurve,
    autocorrelation_decay,
    covariance,
    deviation_correlogram,
    evaluate,
    fit_segments,
    pooled_correlation,
    segment_count
@everywhere using LinearAlgebra: Symmetric, cholesky, diag, dot
@everywhere using StableRNGs: StableRNG
@everywhere using Statistics: mean

@everywhere begin
    const MASSES = collect(126:174)
    const ρ_TRUE = 0.8
    const σ_QUOTED = 0.004
    const MAX_SEGMENTS = 6
    const SEED = 20261127

    # The curve every simulated dataset follows: three segments, joined at 130 and 150.
    function true_ratio(A)
        A ≤ 130 && return 0.5 - 0.08 * (A - 126)
        A ≤ 150 && return 0.18 + 0.021 * (A - 130)
        return 0.60 + 0.25 * (A - 150) / 24
    end

    # The weights of the mean that stands for ⟨R_T⟩: a yield peak at A_H = 140, 6 u wide.
    const YIELD_WEIGHTS = let weight = [exp(-((A - 140) / 6)^2 / 2) for A in MASSES]
        weight ./ sum(weight)
    end

    # One dataset and the measured fraction of its points. An interpolated one is measured on a
    # grid of 1.5 u and written at every mass number, each row counting for the share of a
    # measured point that the retrieval of such a dataset records.
    function simulated_dataset(rng, label, interpolated::Bool, model)
        grid = interpolated ? collect(126:1.5:174) : Float64.(MASSES)
        correlation = Symmetric([ρ_TRUE^abs(a - b) for a in grid, b in grid])
        error =
            model.correlated .* (cholesky(correlation).L * randn(rng, length(grid))) .+
            model.uncorrelated .* randn(rng, length(grid))
        measured = true_ratio.(grid) .+ error .+ σ_QUOTED .* randn(rng, length(grid))
        σ = fill(σ_QUOTED, length(MASSES))
        interpolated || return (RatioCurve(MASSES, measured, σ, label), ones(length(MASSES)))
        written = map(MASSES) do A
            below = searchsortedlast(grid, A)
            below == length(grid) && return measured[end]
            t = (A - grid[below]) / (grid[below + 1] - grid[below])
            return (1 - t) * measured[below] + t * measured[below + 1]
        end
        fraction = fill(length(grid) / length(MASSES), length(MASSES))
        return (RatioCurve(MASSES, written, σ, label), fraction)
    end

    # The trend of a combination, fitted as a run fits it.
    function fitted_trend(pooled, correlation)
        return fit_segments(
            pooled.curve.A_H,
            pooled.curve.ratio,
            pooled.σ_measurement;
            max_segments = MAX_SEGMENTS,
            measured = pooled.measured,
            min_points_per_segment = 4,
            bounds = (0.0, 1.0),
            correlation = correlation,
        )
    end

    # One replicate: the error of the fitted curve and of its weighted mean, and the three
    # uncertainties of each.
    function replicate(seed, pool, model)
        rng = StableRNG(seed)
        sets = [
            simulated_dataset(rng, "set $(i)", pool.interpolated, model) for i in 1:(pool.k)
        ]
        curves = first.(sets)
        fractions = last.(sets)
        pooled = FissionTemperatureRatio._consensus(curves, "systematic_trend", fractions)
        correlogram = deviation_correlogram(curves, pooled.curve; weights = mean.(fractions))
        ρ = autocorrelation_decay(correlogram; lags = 4)
        fit = fitted_trend(
            pooled, pooled_correlation(pooled.curve.A_H, pooled.members, pooled.shares, ρ)
        )
        r_ν = evaluate(fit, MASSES)
        Σ = Matrix(covariance(fit, MASSES))

        # A refit takes its points as independent: only its values are used.
        refits = map(1:(pool.k)) do left_out
            others = setdiff(1:(pool.k), left_out)
            subset = FissionTemperatureRatio._consensus(
                curves[others], "systematic_trend", fractions[others]
            )
            return evaluate(fitted_trend(subset, nothing), MASSES)
        end
        Σ_J = FissionTemperatureRatio._jackknife_covariance(MASSES, refits)

        w = YIELD_WEIGHTS
        σ = (; propagated = sqrt.(max.(diag(Σ), 0.0)), jackknife = sqrt.(diag(Σ_J)))
        σ_mean = (; propagated = sqrt(max(dot(w, Σ, w), 0.0)), jackknife = sqrt(dot(w, Σ_J, w)))
        return (;
            mean = dot(w, r_ν.ratio),
            error = r_ν.ratio .- true_ratio.(MASSES),
            σ = (; σ..., larger = max.(σ.propagated, σ.jackknife)),
            σ_mean = (; σ_mean..., larger = max(σ_mean.propagated, σ_mean.jackknife)),
            segments = segment_count(fit),
            ρ,
        )
    end
end

# The standard deviations of the correlated and of the uncorrelated part of a dataset's error.
function error_model(noise::AbstractString)
    noise == "ar1" && return (; correlated = 0.03, uncorrelated = 0.0)
    noise == "nugget" && return (; correlated = 0.026, uncorrelated = 0.015)
    return throw(ArgumentError("the error model is ar1 or nugget, got $(repr(noise))"))
end

# A pool as the command line names it: the number of datasets and whether they are interpolated.
function parse_pool(pool::AbstractString)
    matched = match(r"^(\d+)([wi])$", pool)
    matched === nothing && throw(ArgumentError("a pool is <k>w or <k>i, got $(repr(pool))"))
    k = parse(Int, matched[1])
    k ≥ 3 || throw(ArgumentError("a pool takes at least three datasets, got $(k)"))
    return (; k, interpolated = matched[2] == "i")
end

function calibrate(name::AbstractString, replicates::Integer, model)
    pool = parse_pool(name)
    rows = pmap(i -> replicate(SEED + 1000 * pool.k + i, pool, model), 1:replicates)
    scatter = std([row.mean for row in rows])
    # The relative standard error of a standard deviation of that many normal values.
    resolution = 1 / sqrt(2 * (replicates - 1))
    for (position, estimate) in enumerate((:propagated, :jackknife, :larger))
        ratio = scatter / sqrt(mean(row.σ_mean[estimate]^2 for row in rows))
        z = reduce(vcat, [row.error ./ row.σ[estimate] for row in rows])
        @printf(
            "%-5s %5s %6s  %-10s  %5.3f ± %5.3f  %6.3f  %6.3f  %6.3f\n",
            position == 1 ? name : "",
            position == 1 ? @sprintf("%.2f", mean(row.segments for row in rows)) : "",
            position == 1 ? @sprintf("%.3f", mean(row.ρ for row in rows)) : "",
            estimate,
            ratio,
            ratio * resolution,
            sqrt(mean(abs2, z)),
            mean(abs.(z) .< 1),
            mean(abs.(z) .< 2),
        )
    end
    return flush(stdout)
end

function main(arguments::Vector{String})
    if length(arguments) < 4 || tryparse(Int, arguments[1]) === nothing
        println(USAGE)
        return 1
    end
    model = error_model(arguments[2])
    replicates = parse(Int, arguments[3])
    replicates ≥ 2 || throw(ArgumentError("at least two replicates, got $(replicates)"))
    println("error model $(arguments[2]): correlated part $(model.correlated) with ρ = \
             $(ρ_TRUE), uncorrelated part $(model.uncorrelated), quoted σ = $(σ_QUOTED); \
             $(replicates) replicates")
    println("mean: standard deviation over the replicates / rms of its σ, 1 calibrated; \
             z = (r̂ − r)/σ at every mass number")
    @printf(
        "%-5s %5s %6s  %-10s  %-13s  %6s  %6s  %6s\n",
        "pool",
        "segs",
        "ρ",
        "σ",
        "mean",
        "rms z",
        "|z|<1",
        "|z|<2",
    )
    foreach(pool -> calibrate(pool, replicates, model), arguments[4:end])
    return 0
end

exit(main(ARGS))
