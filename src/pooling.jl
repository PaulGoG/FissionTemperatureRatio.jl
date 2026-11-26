# The systematic trend: which datasets are pooled, how the datasets of one experiment are
# combined, the fit to the combination and the refits with one dataset left out.

"""
    pooled_datasets(result) -> Vector{String}

The labels of the datasets the systematic trend is combined from, in the order read: those that
provide at least one complete fragment pair within the fragmentation range, that the
configuration does not exclude, and that no republication among them supersedes. A dataset that
forms no pair has nothing to pool, whether or not an exclusion names it.
"""
function pooled_datasets(result::ExtractionResult)
    candidates = _candidates(result.configuration, result.datasets, result.r_ν)
    return [
        result.datasets[i].label for
        i in candidates if !haskey(result.superseded, result.datasets[i].label)
    ]
end

# The indices of the datasets that form a pair in the range and that the configuration does not
# exclude.
function _candidates(
    configuration::Configuration, datasets::Vector{Multiplicity}, r_ν::Vector{RatioCurve}
)
    excluded = _excluded_labels(configuration.excluded_datasets, datasets)
    return [
        i for
        i in eachindex(datasets) if !isempty(r_ν[i]) && !haskey(excluded, datasets[i].label)
    ]
end

"""
    superseded_datasets(result) -> Dict{String,String}

The datasets of a run that a republication supersedes, by label, each with the label of the
dataset that supersedes it. The retrieval marks the earlier publication of one result with the
qualifier `superseded`; it is read, fitted and written, and enters the pool through its successor
alone. Where the successor is excluded or forms no fragment pair, neither is pooled: a withdrawn
result does not stand in for the one that replaced it.
"""
superseded_datasets(result::ExtractionResult) = result.superseded

"""
    deviation_correlogram(result) -> Vector{Union{Missing,Float64}}

The correlogram of the pooled datasets' deviations from the combined curve of a run, lags one to
[`CORRELOGRAM_LAGS`](@ref).
"""
deviation_correlogram(result::ExtractionResult) = result.deviation_correlogram

"""
    deviation_autocorrelation(result) -> Union{Float64,Missing}

The coefficient `ρ` with which the error of a pooled dataset is taken to be correlated along the
mass axis, `ρ^|A − A′|` between two of its mass numbers: [`autocorrelation_decay`](@ref) of
[`deviation_correlogram`](@ref) over the configured number of lags,
`segments.autocorrelation_lags`.

A dataset departs from the others by a slow drift along the mass axis rather than point by point;
with independent points the uncertainty of a smooth curve through the combined values would be
low by about `√((1 + ρ)/(1 − ρ))`. The covariance of the trend, and with it the uncertainty of
every tabulated `R_T` of the trend and of its total average, is formed with the correlation
matrix [`pooled_correlation`](@ref) builds from `ρ` and from the datasets combined at each mass
number. The curve of a single dataset is fitted to its points as independent ones.

`missing` where no lag could be estimated, the trend then being fitted with independent points.
"""
deviation_autocorrelation(result::ExtractionResult) = result.autocorrelation

"""
    leave_one_out_spread(result, mass_yield) -> Union{NamedTuple,Nothing}
    leave_one_out_spread(refits, mass_yield) -> Union{NamedTuple,Nothing}

The spread of the systematic trend's `⟨R_T⟩` over the yield distribution labelled `mass_yield`
across the [`LeaveOneOut`](@ref) refits of a run: `(; min, max, uncertainty, refits)`, the least
and the greatest of the values `θᵢ`, `i = 1, …, k`, of the refits that hold one, the delete-one
jackknife standard error over the pooled datasets or experiments,

```
σ = [(k − 1)/k Σᵢ (θᵢ − θ̄)²]^(1/2),
```

and `k`. This is the uncertainty of the trend's `⟨R_T⟩` that comes from the datasets differing
from one another by more than their errors along the mass axis, which the covariance of the trend
does not contain. `nothing` where fewer than two refits hold a value.

A refit that gave no curve holds no value and is not among the `k`. With two refits the standard
error is half their difference, `|θ₁ − θ₂|/2`: each refit then rests on a single dataset or
experiment, and the number says how far the two lie apart, not how well their trend is known.
"""
function leave_one_out_spread(result::ExtractionResult, mass_yield::AbstractString)
    return leave_one_out_spread(result.leave_one_out, mass_yield)
end

function leave_one_out_spread(refits::Vector{LeaveOneOut}, mass_yield::AbstractString)
    θ = Float64[
        entry.total_average_R_T[mass_yield] for
        entry in refits if haskey(entry.total_average_R_T, mass_yield)
    ]
    k = length(θ)
    k < 2 && return nothing
    θ̄ = mean(θ)
    return (;
        min = minimum(θ),
        max = maximum(θ),
        uncertainty = sqrt((k - 1) / k * sum(abs2, θ .- θ̄)),
        refits = k,
    )
end

# The trend refitted with each pooled dataset left out in turn, the datasets of one experiment
# together: how far ⟨R_T⟩, the number of segments and the breakpoints rest on one of them. Only the
# values of a refit are used, so it takes the points as independent. The curves of the refits
# are returned beside their records; the jackknife of the trend is formed from them.
function _leave_one_out(
    units,
    datasets::Vector{Multiplicity},
    relation,
    configuration::Configuration,
    averaged::Vector{MassYield},
)
    length(units.curves) ≥ 2 || return (; entries = LeaveOneOut[], curves = ExtractedCurve[])
    refits = LeaveOneOut[]
    curves = ExtractedCurve[]
    settings = configuration.segments
    @info "leave-one-out refits of the systematic trend" refits = length(units.curves)
    for j in eachindex(units.curves)
        others = setdiff(eachindex(units.curves), j)
        refitted = _fit_trend(
            units.curves[others],
            units.weights[others],
            relation,
            settings,
            configuration.system.compound.A;
            estimate_autocorrelation = false,
            quiet = true,
        )
        left_out = [datasets[i].label for i in units.members[j]]
        refit = refitted.trend
        if !(refit isa ExtractedCurve)
            push!(
                refits,
                LeaveOneOut(
                    left_out, "no fit: $(refit.reason)", 0, Int[], NaN, Dict{String,Float64}()
                ),
            )
            continue
        end
        outcome = if refitted.windows == settings.required_windows
            "segmented curve"
        else
            "segmented curve, fitted without the required windows"
        end
        averages = Dict{String,Float64}(
            distribution.label => first(total_average(refit, distribution)) for
            distribution in averaged if _overlaps(refit.R_T, distribution)
        )
        push!(curves, refit)
        push!(
            refits,
            LeaveOneOut(
                left_out,
                outcome,
                segment_count(refit.fit),
                copy(refit.fit.breakpoints),
                refit.fit.wrss / refit.fit.dof,
                averages,
            ),
        )
    end
    # What the jackknife over the refits can and cannot say here.
    failed = [join(entry.datasets, " + ") for entry in refits if entry.segments == 0]
    isempty(failed) ||
        @warn "a leave-one-out refit gave no curve; the jackknife uncertainty of the trend is \
               formed over the other refits" left_out = failed
    length(units.curves) == 2 &&
        @warn "two pooled experiments: each leave-one-out refit rests on one of them, and the \
               jackknife uncertainty of the trend is half the difference of the two"
    return (; entries = refits, curves)
end

# The experiments among the datasets `candidates` indexes: the connected components of the
# relation their retrieval records state under `correlated_with`, each a list of positions in
# `candidates`, in the order of their first members.
function _experiments(datasets::Vector{Multiplicity}, candidates::Vector{Int})
    accessions = [_accession(datasets[i]) for i in candidates]
    partners = [_correlated_with(datasets[i].record) for i in candidates]
    experiment = collect(eachindex(candidates))
    for a in eachindex(candidates), b in (a + 1):length(candidates)
        named =
            (!isempty(accessions[b]) && accessions[b] in partners[a]) ||
            (!isempty(accessions[a]) && accessions[a] in partners[b])
        named || continue
        merged, kept = experiment[b], experiment[a]
        replace!(experiment, merged => kept)
    end
    return [findall(==(first_position), experiment) for first_position in unique(experiment)]
end

# How the datasets of one experiment are related, where every one of their records states it and
# they agree; empty where a member states none or two members differ. A relation is a statement
# about the experiment, and one its records do not share is not acted on.
function _relation(datasets::Vector{Multiplicity}, members::Vector{Int})
    stated = [_correlation_relation(datasets[i].record) for i in members]
    return (all(!isempty, stated) && allequal(stated)) ? first(stated) : ""
end

# The datasets among `candidates` that a republication supersedes, each with the dataset that
# supersedes it: in an experiment whose records give the relation `republication`, the members
# carrying the qualifier `superseded`, where exactly one member does not. Such a dataset is the
# same result published earlier and enters no pool: not beside its successor, and not in its
# place where the successor is excluded or forms no fragment pair, since a withdrawn result does
# not stand in for the one that replaced it. The experiments are those of every dataset read, so
# that a successor outside the pool is still seen.
function _superseded(datasets::Vector{Multiplicity}, candidates::Vector{Int})
    superseded = Dict{Int,Int}()
    for members in _experiments(datasets, collect(eachindex(datasets)))
        length(members) > 1 || continue
        _relation(datasets, members) == "republication" || continue
        withdrawn = [i for i in members if "superseded" in _qualifiers(datasets[i].record)]
        current = setdiff(members, withdrawn)
        length(current) == 1 || continue
        for i in withdrawn
            i in candidates && (superseded[i] = only(current))
        end
    end
    return superseded
end

# Two or more reductions of one set of events, combined mass number by mass number. Their
# statistical errors are common, so nothing is gained by averaging them: the value is their mean,
# and its uncertainty the largest a member quotes with half the difference between the members
# added in quadrature, which carries the uncertainty of the reduction. A mass number one member
# alone holds takes that member's value and uncertainty; where no member quotes an uncertainty
# the point quotes none. Returned with the measured fraction of each point, which divides the
# quoted part alone: with f̄ the mean fraction of the members there, (σ² + d²)/(σ²/f̄ + d²),
# σ the largest uncertainty quoted and d the half difference. Interpolation spreads a
# statistical error over the rows written; the difference between two reductions is not one.
function _alternative_analyses(
    curves::Vector{RatioCurve}, fractions::AbstractVector{<:Real}, label::AbstractString
)
    masses = sort!(unique!(reduce(vcat, (curve.A_H for curve in curves); init = Int[])))
    ratio = Float64[]
    σ = Union{Missing,Float64}[]
    measured = Float64[]
    for mass in masses
        values = Float64[]
        quoted = Float64[]
        held = Float64[]
        for (curve, fraction) in zip(curves, fractions)
            index = findfirst(==(mass), curve.A_H)
            index === nothing && continue
            push!(values, curve.ratio[index])
            ismissing(curve.σ[index]) || push!(quoted, curve.σ[index])
            push!(held, fraction)
        end
        half_difference = (maximum(values) - minimum(values)) / 2
        push!(ratio, mean(values))
        push!(σ, isempty(quoted) ? missing : hypot(maximum(quoted), half_difference))
        f̄ = mean(held)
        if isempty(quoted)
            push!(measured, f̄)
        else
            s = maximum(quoted)
            d = half_difference
            total = s^2 + d^2
            push!(measured, total == 0 ? f̄ : total / (s^2 / f̄ + d^2))
        end
    end
    return RatioCurve(masses, ratio, σ, String(label)), measured
end

# What a trend is combined from: each admitted dataset, or, for the datasets of one experiment —
# those whose retrieval records name one another under `correlated_with` — one curve that stands
# for the experiment, entering the pool at the uncertainty of one measurement with the measured
# fraction of each of its points. Runs or analyses on one apparatus share their errors; pooled
# side by side they would count as independent measurements in the between-dataset variance and
# weigh as several. How the curve is formed follows the relation the records state: alternative
# analyses of the same events by `_alternative_analyses`; repeated runs, parts of one spectrum
# and members whose records state no relation, or different ones, as a pool is combined. Of a
# republication the superseding dataset alone is admitted, `_superseded`. For each curve,
# `weights` holds the measured fraction of every point and `members` the indices of its datasets.
function _pool_units(
    datasets::Vector{Multiplicity}, r_ν::Vector{RatioCurve}, admitted::Vector{Int}
)
    curves = RatioCurve[]
    weights = Vector{Float64}[]
    members = Vector{Int}[]
    for positions in _experiments(datasets, admitted)
        indices = admitted[positions]
        push!(members, indices)
        fractions = [_measured_fraction(datasets[i].record) for i in indices]
        if length(indices) == 1
            curve = r_ν[only(indices)]
            push!(curves, curve)
            push!(weights, fill(only(fractions), length(curve)))
            continue
        end
        label = join((datasets[i].label for i in indices), " + ")
        if _relation(datasets, indices) == "alternative_analysis"
            curve, measured = _alternative_analyses(r_ν[indices], fractions, label)
            push!(curves, curve)
            push!(weights, measured)
        else
            combination = _consensus(
                r_ν[indices],
                label,
                [fill(f, length(r_ν[i])) for (i, f) in zip(indices, fractions)],
            )
            push!(
                curves,
                RatioCurve(
                    combination.curve.A_H,
                    combination.curve.ratio,
                    combination.σ_measurement,
                    label,
                ),
            )
            push!(weights, combination.measured)
        end
    end
    return (; curves, weights, members)
end

# The systematic trend of a pool of datasets: their combination, the correlogram of their
# deviations from it, the autocorrelation fitted to that, and the segmented curve, fitted with the
# required windows and, where no curve satisfies them, without. `combined` is the combination at
# the uncertainty of one measurement, the curve the trend is fitted to; the trend is the curve, or
# the reason there is none; `windows` are those it was fitted with, and `autocorrelation` is
# `missing` where no lag of the correlogram could be estimated. Without `estimate_autocorrelation`
# the points are taken as independent, which suffices where only the values of the fit are used.
function _fit_trend(
    curves::Vector{RatioCurve},
    weights::Vector{Vector{Float64}},
    relation::Tuple{RatioAveraging,LevelDensityModel,FragmentationDomain},
    settings::SegmentSettings,
    A₀::Integer;
    estimate_autocorrelation::Bool = true,
    quiet::Bool = false,
)
    pooled = _consensus(curves, SYSTEMATIC_TREND_LABEL, weights)
    # The trend is fitted to the combined values at the uncertainty of one measurement, with the
    # measured fraction beside it, as a dataset is: the standard error of the combination carries
    # the fraction already, and fitting to it would count the fraction twice.
    combined = RatioCurve(
        pooled.curve.A_H, pooled.curve.ratio, pooled.σ_measurement, pooled.curve.label
    )
    correlogram = if estimate_autocorrelation
        # A curve that is itself a combination enters the correlogram with its mean fraction.
        deviation_correlogram(curves, pooled.curve; weights = mean.(weights))
    else
        Vector{Union{Missing,Float64}}(missing, CORRELOGRAM_LAGS)
    end
    lags = min(settings.autocorrelation_lags, length(correlogram))
    estimated = !all(ismissing, view(correlogram, 1:lags))
    ρ = estimated ? autocorrelation_decay(correlogram; lags = lags) : 0.0
    # The errors of the combined points: correlated along the mass axis within a dataset,
    # independent between datasets.
    correlation = if ρ > 0
        pooled_correlation(pooled.curve.A_H, pooled.members, pooled.shares, ρ)
    else
        nothing
    end
    fit(windows) = _segment(
        combined,
        relation,
        settings,
        SYSTEMATIC_TREND_LABEL,
        windows,
        A₀;
        measured = pooled.measured,
        correlation = correlation,
        quiet = quiet,
    )
    windows = settings.required_windows
    trend = fit(windows)
    if trend isa NoFit && !isempty(windows)
        quiet ||
            @warn "no segmented curve satisfies the required windows; the systematic-trend \
                   curve was fitted without them" windows reason = trend.reason
        windows = UnitRange{Int}[]
        trend = fit(windows)
    end
    autocorrelation = estimated ? ρ : missing
    return (; pooled, combined, correlogram, autocorrelation, trend, windows)
end
