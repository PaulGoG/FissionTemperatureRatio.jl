# Calibration of the uncertainty of the systematic trend on the pool of a configuration:
#
#     julia scripts/trend_uncertainty_calibration.jl <workers> <configuration.toml> <replicates> [<seed>]
#
# The pooled datasets of the configuration, the datasets of one experiment as one, keep their mass
# numbers, their quoted uncertainties and their measured fractions, and take simulated values
# about a known curve, the trend their own values give:
#
#     y_i(A) = r(A) + s_i e_i(A) + σ_i(A) ε_i(A).
#
# r is the fitted multiplicity ratio of that trend; σ_i is the uncertainty dataset i quotes, none
# where it quotes none; ε_i is independent and standard normal; e_i has unit variance and the
# correlation ρ^|A − A′| along the mass axis, ρ being the coefficient the run estimates; and
# s_i² = max(⟨(y_i − r)²⟩ − ⟨σ_i²⟩, 0) is what the residuals of dataset i about the trend hold
# beyond its quoted uncertainty. A replicate is then treated as a run treats the data: the pool
# is combined, ρ is fitted to the correlogram of the deviations, the trend is fitted with the pin,
# the required windows and the selection of the configuration, and it is refitted with half of the
# pooled datasets left out, over the subsets `segments.jackknife_subsets` and
# `segments.jackknife_seed` give. Three uncertainties of R_T follow: the propagated one, the
# jackknife over the refits, and the larger of the two, which is what a run reports.
#
# Scored over the replicates, against the R_T of the known curve, at every mass number but a pinned
# one: the root mean square of z = (R̂_T − R_T)/σ and the fractions of |z| below 1 and below 2,
# which are 0.683 and 0.954 for a normal error with a known σ; and the 5 %, 50 % and 95 % quantiles
# of σ over the root-mean-square error at its mass number, which say how far a single reported
# uncertainty lies from the error it stands for. For ⟨R_T⟩ over the primary yield distribution, or
# the unweighted mean over the range without one: the standard deviation over the replicates
# divided by the root mean square of the uncertainty, 1 for a calibrated uncertainty and above 1
# for one that is too small, with the resolution 1/√(2(n − 1)) the number of replicates gives it,
# and the 5 %, 50 % and 95 % quantiles of the uncertainty over that standard deviation.
#
# The replicates are shared among `<workers>` processes, or run in this one for 0. A replicate
# costs one trend fit and one per refit: about ten seconds for four pooled experiments and a few
# minutes for thirteen. The seed of each replicate is drawn from one stream seeded by `<seed>`,
# 20261127 where none is given, so the same arguments print the same table whatever the number of
# workers. The first lines printed record the arguments and the design. The pooling, the fit and
# the jackknife are the package's own internal functions, so that what is scored is what a run
# computes. Nothing is written.

include(joinpath(@__DIR__, "activate.jl"))

using Distributed: @everywhere, addprocs, pmap
using Printf: @printf
using Statistics: mean, std

const USAGE = "usage: julia scripts/trend_uncertainty_calibration.jl <workers> \
               <configuration.toml> <replicates> [<seed>]"
const ROOT = dirname(@__DIR__)

# The arguments, or `nothing` where they are not ones the script can run on.
function parsed(arguments::Vector{String})
    length(arguments) in (3, 4) || return nothing
    workers = tryparse(Int, arguments[1])
    replicates = tryparse(Int, arguments[3])
    seed = length(arguments) == 4 ? tryparse(Int, arguments[4]) : 20261127
    (workers === nothing || replicates === nothing || seed === nothing) && return nothing
    (workers ≥ 0 && replicates ≥ 2 && seed ≥ 0) || return nothing
    # A path is taken as given when it exists from the working directory, and otherwise relative
    # to the repository root, as `scripts/run.jl` takes it.
    path = isfile(arguments[2]) ? abspath(arguments[2]) : joinpath(ROOT, arguments[2])
    isfile(path) || return nothing
    return (; workers, path, replicates, seed)
end

const SETTINGS = parsed(ARGS)
if SETTINGS === nothing
    println(USAGE)
    exit(1)
end

# The workers run in this environment, which the activation above has resolved already.
SETTINGS.workers > 0 && addprocs(SETTINGS.workers; exeflags = `--project=$(@__DIR__)`)

@everywhere using FissionFragmentsDomain: mass_yield
@everywhere using FissionTemperatureRatio:
    FissionTemperatureRatio, ExtractedCurve, RatioCurve, evaluate, load_configuration
@everywhere using LinearAlgebra: Symmetric, cholesky, diag, dot
@everywhere using Logging: NullLogger, with_logger
@everywhere using StableRNGs: StableRNG
@everywhere using Statistics: mean

@everywhere begin
    const FTR = FissionTemperatureRatio

    # The pool of the configuration and the simulation built on it: the known curve at the mass
    # numbers of every pooled unit, what each quotes, what its residuals hold beyond that, and the
    # factor of its correlation matrix.
    function staged(path::AbstractString, data::AbstractString)
        configuration = load_configuration(path; data_directory = data)
        return with_logger(NullLogger()) do
            retrievals = FTR._read_retrievals(FTR._retrieval_directories(configuration))
            inversion = FTR._inversion(configuration, retrievals)
            relation = (inversion.averaging, inversion.model, inversion.domain)
            input = FTR._multiplicity_input(configuration, inversion, retrievals)
            pool = FTR._pooled_trend(configuration, relation, input)
            trend = pool.fitted.trend
            trend isa ExtractedCurve ||
                error("the configuration gives no systematic trend: $(trend.reason)")
            units = pool.units
            k = length(units.curves)
            k ≥ 2 || error("one pooled dataset: there is nothing to leave out")
            ρ = coalesce(pool.fitted.autocorrelation, 0.0)
            known = [
                Float64[first(evaluate(trend.fit, A)) for A in c.A_H] for c in units.curves
            ]
            quoted = [Float64[coalesce(σ, 0.0) for σ in c.σ] for c in units.curves]
            unquoted = [
                sqrt(max(mean(abs2, c.ratio .- r) - mean(abs2, q), 0.0)) for
                (c, r, q) in zip(units.curves, known, quoted)
            ]
            factors = [
                cholesky(Symmetric([ρ^abs(a - b) for a in c.A_H, b in c.A_H])).L for
                c in units.curves
            ]
            yields = FTR._mass_yield_input(configuration, inversion, retrievals)
            masses = trend.R_T.A_H
            weight = if isempty(yields.averaged)
                ones(length(masses))
            else
                primary = first(yields.averaged)
                Float64[something(mass_yield(primary, A), (0.0, 0.0))[1] for A in masses]
            end
            d = FTR._jackknife_deleted(k)
            settings = configuration.segments
            chosen = FTR._left_out_subsets(
                k, d, settings.jackknife_subsets, settings.jackknife_seed
            )
            return (;
                configuration,
                relation,
                units,
                k,
                d,
                chosen,
                ρ,
                known,
                quoted,
                unquoted,
                factors,
                masses,
                truth = trend.R_T.ratio,
                pinned = trend.fit.pinned_value !== nothing,
                weight = weight ./ sum(weight),
                weighted = !isempty(yields.averaged),
            )
        end
    end

    # The trend of the units `kept` of `curves`, fitted as a run fits it and, for a refit, as a
    # run fits its refits, with the autocorrelation fitted; `nothing` where there is no trend.
    function fitted(stage, curves, kept; refit::Bool)
        result = with_logger(NullLogger()) do
            FTR._fit_trend(
                curves[kept],
                stage.units.weights[kept],
                stage.relation,
                stage.configuration.segments,
                stage.configuration.system.compound.A;
                estimate_autocorrelation = !refit,
                quiet = true,
            )
        end
        result.trend isa ExtractedCurve || return nothing
        return (; trend = result.trend, ρ = coalesce(result.autocorrelation, 0.0))
    end

    # One replicate: the errors of R_T and of its weighted mean, with the three uncertainties of
    # each; `nothing` where the replicate gives no trend over the mass numbers of the known curve
    # or fewer than two refits.
    function replicate(stage, seed::UInt64)
        rng = StableRNG(seed)
        curves = map(eachindex(stage.units.curves)) do i
            curve = stage.units.curves[i]
            n = length(curve.A_H)
            values =
                stage.known[i] .+ stage.unquoted[i] .* (stage.factors[i] * randn(rng, n)) .+
                stage.quoted[i] .* randn(rng, n)
            return RatioCurve(curve.A_H, values, curve.σ, curve.label)
        end
        everything = collect(1:(stage.k))
        whole = fitted(stage, curves, everything; refit = false)
        (whole === nothing || whole.trend.R_T.A_H != stage.masses) && return nothing
        trend = whole.trend
        refits = ExtractedCurve[]
        for subset in stage.chosen.subsets
            refit = fitted(stage, curves, setdiff(everything, subset); refit = true)
            refit === nothing || push!(refits, refit.trend)
        end
        length(refits) ≥ 2 || return nothing
        reported = ExtractedCurve(trend, refits; units = stage.k, deleted = stage.d)
        Σ = trend.R_T_covariance
        Σ_J = reported.R_T_jackknife_covariance
        w = stage.weight
        propagated = sqrt(max(dot(w, Σ, w), 0.0))
        jackknife = sqrt(max(dot(w, Σ_J, w), 0.0))
        return (;
            error = trend.R_T.ratio .- stage.truth,
            σ = (Float64.(trend.R_T.σ), sqrt.(max.(diag(Σ_J), 0.0)), Float64.(reported.R_T.σ)),
            average = dot(w, trend.R_T.ratio) - dot(w, stage.truth),
            σ_average = (propagated, jackknife, max(propagated, jackknife)),
            segments = FTR.segment_count(trend.fit),
            refits = length(refits),
            ρ = whole.ρ,
        )
    end
end

const DATA = joinpath(ROOT, "data")
@everywhere const STAGE = staged($(SETTINGS.path), $DATA)
@everywhere job(seed::UInt64) = replicate(STAGE, seed)

const NAMES = ("propagated", "jackknife", "reported")

function quantile_of(values::Vector{Float64}, p::Real)
    return sort(values)[clamp(round(Int, p * length(values)), 1, length(values))]
end

function main()
    stage = STAGE
    seeds = rand(StableRNG(SETTINGS.seed), UInt64, SETTINGS.replicates)
    elapsed = @elapsed outcomes = pmap(job, seeds)
    rows = filter(!isnothing, outcomes)
    length(rows) ≥ 2 || error("fewer than two replicates gave a trend with its refits")

    @printf(
        "%s: %d replicates, seed %d, %d workers, %.1f min\n",
        relpath(SETTINGS.path, ROOT),
        SETTINGS.replicates,
        SETTINGS.seed,
        SETTINGS.workers,
        elapsed / 60
    )
    @printf(
        "  %d pooled datasets or experiments, %d left out of each of %d refits (%s); ρ = %.3f\n",
        stage.k,
        stage.d,
        length(stage.chosen.subsets),
        stage.chosen.complete ? "every subset" : "drawn of $(binomial(big(stage.k), stage.d))",
        stage.ρ
    )
    println("  error beyond the quoted uncertainty, by pooled dataset or experiment:")
    for (curve, s) in zip(stage.units.curves, stage.unquoted)
        @printf("    %-50s %.4f\n", curve.label, s)
    end
    @printf(
        "  %d replicates gave a trend over the range with its refits; on average %.2f \
         segments selected, ρ fitted as %.3f, %.1f of %d refits with a curve\n",
        length(rows),
        mean(row.segments for row in rows),
        mean(row.ρ for row in rows),
        mean(row.refits for row in rows),
        length(stage.chosen.subsets)
    )

    scored = (stage.pinned ? 2 : 1):length(stage.masses)
    rmse = [sqrt(mean(abs2, [row.error[j] for row in rows])) for j in eachindex(stage.masses)]
    scatter = std([row.average for row in rows])
    resolution = 1 / sqrt(2 * (length(rows) - 1))
    @printf(
        "  %s: scatter over the replicates %.4f\n",
        stage.weighted ? "⟨R_T⟩ over the primary yield" : "the unweighted mean of R_T",
        scatter
    )
    @printf(
        "  %-11s %7s %9s %9s   %s   %s   %s\n",
        "uncertainty",
        "rms z",
        "P(|z|<1)",
        "P(|z|<2)",
        "σ / rms error: 5 %   50 %   95 %",
        "mean: scatter / rms σ",
        "σ / scatter: 5 %   50 %   95 %"
    )
    for (position, name) in enumerate(NAMES)
        z = Float64[]
        relative = Float64[]
        # An uncertainty of zero beside an error, as the jackknife gives at a mass number fewer
        # than two refits hold, covers nothing: it counts in the fractions and not in rms z.
        zeros_found = 0
        for row in rows, j in scored
            σ = row.σ[position][j]
            σ > 0 ? push!(z, row.error[j] / σ) : (zeros_found += 1)
            push!(relative, σ / rmse[j])
        end
        total = length(z) + zeros_found
        calibration = scatter / sqrt(mean(abs2, [row.σ_average[position] for row in rows]))
        of_mean = [row.σ_average[position] / scatter for row in rows]
        @printf(
            "  %-11s %7.3f %9.3f %9.3f                  %5.2f  %5.2f  %5.2f        %5.3f ± %5.3f                %5.2f  %5.2f  %5.2f%s\n",
            name,
            sqrt(mean(abs2, z)),
            count(<(1), abs.(z)) / total,
            count(<(2), abs.(z)) / total,
            quantile_of(relative, 0.05),
            quantile_of(relative, 0.5),
            quantile_of(relative, 0.95),
            calibration,
            calibration * resolution,
            quantile_of(of_mean, 0.05),
            quantile_of(of_mean, 0.5),
            quantile_of(of_mean, 0.95),
            zeros_found == 0 ? "" : "   ($(zeros_found) of $(total) uncertainties are zero)"
        )
    end
    return 0
end

exit(main())
