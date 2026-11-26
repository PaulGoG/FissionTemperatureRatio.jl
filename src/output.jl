# Tabulated output: the tables and the manifest of a run, and the names of its files.

"""
    curve_accessions(result) -> Dict{String,String}

The EXFOR accession of every dataset of a run, keyed by its label: the identifier the retrieval
record states, or the one the file name leads with where no record lists the file, and empty for a
tabulation that carries neither. The systematic trend is fitted to several datasets and has
none; its entry is empty.

The label, `Author year`, is the display name and the key a consuming code selects a curve by; it
carries the accession only where two datasets would otherwise share it. The accession is what
identifies the measurement in the archive, so it is written into the manifest entry of every
dataset curve and into the name of every file of that dataset.
"""
function curve_accessions(result::ExtractionResult)
    accessions = Dict(data.label => _accession(data) for data in result.datasets)
    accessions[SYSTEMATIC_TREND_LABEL] = ""
    return accessions
end

# The token the tables and figures of each curve are named by, keyed by label. A dataset's is the
# stem of its input file, `<accession>_<Author>_<year>` for a retrieved one, so an output traces to
# its archive entry and to its input file by name alone; the trend's is its label.
function _file_tokens(result::ExtractionResult)
    tokens = Dict(
        data.label => _file_token(first(splitext(basename(data.source)))) for
        data in result.datasets
    )
    tokens[SYSTEMATIC_TREND_LABEL] = _file_token(SYSTEMATIC_TREND_LABEL)
    allunique(values(tokens)) || throw(
        ArgumentError(
            "two datasets would be written under one file name: $(sort!(collect(values(tokens))))",
        ),
    )
    return tokens
end

# The run record a consuming code reads: the system, the quantity tabulated, the domain the curves
# were extracted on, and for every segmented curve its label, its kind, the two files it is
# tabulated in and, for a dataset curve, its archive accession, written by FissionFragmentsDomain's
# writer and nothing besides. What else a run knows about each curve is in the segmented-curves
# table and in the run metadata.
function _write_manifest(
    result::ExtractionResult,
    directory::AbstractString,
    identifier::AbstractString,
    written::Dict{String,String},
)
    record = system_record(result.configuration.system)
    system = ManifestSystem(
        record["label"],
        record["notation"],
        record["target_A"],
        record["target_Z"],
        record["channel"],
        record["reaction"],
        record["incident_energy_MeV"],
        record["compound_A"],
        record["compound_Z"],
    )
    accessions = curve_accessions(result)
    curves = [
        ManifestCurve(
            curve.label,
            _kind(curve),
            basename(written["R_T_segmented/$(curve.label)"]),
            basename(written["r_nu_pivots/$(curve.label)"]);
            accession = accessions[curve.label],
        ) for curve in result.segmented_curves
    ]
    path = joinpath(directory, "manifest_$(identifier).toml")
    manifest = TemperatureRatioManifest(
        system,
        MANIFEST_ORDINATE,
        copy(MANIFEST_ABSCISSA),
        # Stated so that a reader knows what the columns hold, not so that a reader looks them up
        # by name: the files are read by column position.
        ["A_H", "R_T", "R_T_uncertainty"],
        curves,
        path,
        manifest_domain(result),
    )
    return write_temperature_ratio_manifest(path, manifest)
end

# A ratio table: the abscissa, the quantity, then its uncertainty. The quantity names its own
# column, so a reader holding the file knows what is in it without consulting the file name. An
# unquoted uncertainty is an empty field, never a zero, which would denote an exact value.
function _ratio_table(curve::RatioCurve, quantity::AbstractString, significant_digits::Integer)
    return DataFrame(
        :A_H => curve.A_H,
        Symbol(quantity) => round.(curve.ratio; sigdigits = significant_digits),
        Symbol("$(quantity)_uncertainty") => [
            ismissing(σ) ? missing : round(σ; sigdigits = significant_digits) for σ in curve.σ
        ],
    )
end

"""
    write_results(result, directory) -> Dict{String,String}

Write the tabulated ratios, the segment pivots, the segmented-curve summary, the total averages,
the leave-one-out refits of the systematic trend, the dataset diagnostics and the manifest of a
run into `directory`, and return the paths written, keyed by content.

`directory` is created. One that already holds files is refused rather than written into, so a
caller that wants a second run of one configuration beside the first moves the first aside;
`scripts/run.jl` does, numbering it `#1`, `#2`, …, and names the directory
`data/sims/<system>/<run identifier>/`. The provenance record and the figures are not written
here: [`run_metadata`](@ref) supplies what the library knows about a run, and the script adds the
rest, writes `metadata.toml`, and calls [`write_figures`](@ref).

The manifest, `manifest_<run identifier>.toml`, is the contract with a consuming code, which
stages the whole directory and selects the manifest by that prefix; it holds the system, the
domain and, per segmented curve, the label, the kind, the two files and the accession of a dataset
curve, and nothing else. What else a run reports about each curve — segments, pin, span, pairs,
coverage, reduced chi-squared, range mean — is one row per manifest curve, keyed by its label,
in `segmented_curves_<run identifier>.csv`; the total averages are in
`total_average_R_T_<run identifier>.csv`, and the trend refitted with each pooled dataset left
out, where the run made those refits, in `leave_one_out_<run identifier>.csv`.
`dataset_diagnostics.csv` flags, with the reason, every dataset curve that resolves no minimum
over a range covering a required window.

The tables of a dataset are named by the stem of its input file, which for a retrieved dataset is
`<accession>_<Author>_<year>`: `R_T_vs_A_H_segmented_<accession>_<Author>_<year>.csv`. The label
stays `Author year`, and [`curve_accessions`](@ref) gives the accession of each.
"""
function write_results(result::ExtractionResult, directory::AbstractString)
    isdir(directory) &&
        !isempty(readdir(directory)) &&
        throw(
            ArgumentError(
                "$(directory) already holds files; write_results never writes into a directory \
                 that does, move it aside first"
            ),
        )
    identifier = run_identifier(result.configuration)
    tokens = _file_tokens(result)
    # A name longer than a file system admits would fail halfway through writing the run, so
    # every name is checked before anything is written.
    for name in _file_names(identifier, tokens)
        ncodeunits(name) <= MAX_FILE_NAME_BYTES || throw(
            ArgumentError("the file name $(name) is $(ncodeunits(name)) bytes, beyond the \
                 $(MAX_FILE_NAME_BYTES) a file system admits")
        )
    end
    mkpath(directory)
    digits = result.configuration.output.significant_digits
    accessions = curve_accessions(result)
    written = Dict{String,String}()
    table(name, rows) = begin
        path = joinpath(directory, name)
        CSV.write(path, DataFrame(rows))
        path
    end

    _write_ratio_tables!(written, result, directory, tokens, digits)
    written["segmented_curves"] = table(
        "segmented_curves_$(identifier).csv", _segmented_curve_rows(result, accessions, digits)
    )
    averages = _total_average_rows(result, digits)
    isempty(averages) ||
        (written["total_average_R_T"] = table("total_average_R_T_$(identifier).csv", averages))
    isempty(result.leave_one_out) || (
        written["leave_one_out"] = table(
            "leave_one_out_$(identifier).csv",
            _leave_one_out_rows(result, accessions, digits),
        )
    )
    isempty(result.dataset_diagnostics) || (
        written["dataset_diagnostics"] = table(
            "dataset_diagnostics.csv", _diagnostic_rows(result, accessions, digits)
        )
    )
    written["manifest"] = _write_manifest(result, directory, identifier, written)
    @info "results written" directory
    return written
end

# The ratio tables: r_ν and R_T point by point for every dataset, the combination the trend was
# fitted to, and for every segmented curve its pivots, with the uncertainty of the curve at each,
# and its r_ν and R_T at every mass number of its range.
function _write_ratio_tables!(
    written::Dict{String,String},
    result::ExtractionResult,
    directory::AbstractString,
    tokens::AbstractDict,
    digits::Integer,
)
    for (curves, name, quantity) in (
        (result.r_ν, "r_nu_vs_A_H", "r_nu"),
        (result.R_T, "R_T_vs_A_H", "R_T"),
        ([result.consensus_r_ν], "r_nu_vs_A_H_consensus", "r_nu"),
    )
        for curve in curves
            isempty(curve) && continue
            path = joinpath(directory, "$(name)_$(tokens[curve.label]).csv")
            CSV.write(path, _ratio_table(curve, quantity, digits))
            written["$(name)/$(curve.label)"] = path
        end
    end
    for curve in result.segmented_curves
        token = tokens[curve.label]
        points = pivots(curve.fit)
        pivot_table = DataFrame(;
            A_H = [point[1] for point in points],
            r_nu = [round(point[2]; sigdigits = digits) for point in points],
            # The uncertainty of the curve's `r_ν` at each pivot: the propagated one, and for a
            # trend that carries a jackknife the larger of the two. A file that states the
            # segmented description without it invites the reader to assume there is none.
            r_nu_uncertainty = [
                round(curve.r_ν.σ[findfirst(==(point[1]), curve.r_ν.A_H)]; sigdigits = digits)
                for point in points
            ],
        )
        path = joinpath(directory, "r_nu_vs_A_H_pivots_$(token).csv")
        CSV.write(path, pivot_table)
        written["r_nu_pivots/$(curve.label)"] = path
        for (segmented, name, quantity) in (
            (curve.r_ν, "r_nu_vs_A_H_segmented", "r_nu"),
            (curve.R_T, "R_T_vs_A_H_segmented", "R_T"),
        )
            path = joinpath(directory, "$(name)_$(token).csv")
            CSV.write(path, _ratio_table(segmented, quantity, digits))
            written["$(quantity)_segmented/$(curve.label)"] = path
        end
    end
    return written
end

# What a run knows about each segmented curve beyond the manifest: one row per manifest curve,
# keyed by its label. No column is an energy; the ratios are dimensionless.
function _segmented_curve_rows(
    result::ExtractionResult, accessions::AbstractDict, digits::Integer
)
    pooled = Set(pooled_datasets(result))
    autocorrelation = deviation_autocorrelation(result)
    return [
        (
            label = curve.label,
            accession = _text_field(accessions[curve.label]),
            kind = _kind(curve),
            pooled = _is_trend(curve) || curve.label in pooled,
            segments = segment_count(curve.fit),
            pinned_at_symmetric_split = curve.fit.pinned_value !== nothing,
            first_A_H = first(curve.R_T.A_H),
            last_A_H = last(curve.R_T.A_H),
            pairs = curve.fit.points,
            measured_points = round(curve.fit.measured_points; sigdigits = digits),
            coverage = round(curve.coverage; sigdigits = digits),
            reduced_chi_squared = round(curve.fit.wrss / curve.fit.dof; sigdigits = digits),
            chi_squared_over_expectation = round(
                curve.fit.wrss / curve.fit.expected_wrss; sigdigits = digits
            ),
            weights_imputed = curve.fit.weights_imputed,
            range_mean_R_T = round(
                first(result.range_mean_R_T[curve.label]); sigdigits = digits
            ),
            range_mean_R_T_uncertainty = round(
                last(result.range_mean_R_T[curve.label]); sigdigits = digits
            ),
            # The autocorrelation the covariance of the trend was formed with.
            deviation_autocorrelation = if _is_trend(curve) && !ismissing(autocorrelation)
                round(autocorrelation; sigdigits = digits)
            else
                missing
            end,
        ) for curve in result.segmented_curves
    ]
end

# The total average over each yield distribution: one row per segmented curve and distribution,
# which is how the literature tabulates it. The uncertainty of the curve first, the propagated
# one and, for a trend that carries a jackknife, the larger of that and the jackknife; the
# independent-points one, the published approximation, beside it; then the fraction of the
# distribution's yield the curve takes in. For the systematic trend, ⟨R_T⟩ refitted with one and
# two segments more than selected, and, as a diagnostic, the delete-one jackknife of ⟨R_T⟩ over
# the pooled datasets with the extremes of ⟨R_T⟩ with one of them left out.
function _total_average_rows(result::ExtractionResult, digits::Integer)
    rounded(x) = round(x; sigdigits = digits)
    rows = NamedTuple{
        (
            :segmented_curve,
            :mass_yield,
            :R_T,
            :R_T_uncertainty,
            :R_T_uncertainty_independent_points,
            :R_T_uncertainty_leave_one_out,
            :R_T_leave_one_out_min,
            :R_T_leave_one_out_max,
            :yield_fraction,
            :R_T_one_more_segment,
            :R_T_two_more_segments,
        ),
        Tuple{
            String,
            String,
            Float64,
            Float64,
            Float64,
            Union{Missing,Float64},
            Union{Missing,Float64},
            Union{Missing,Float64},
            Float64,
            Union{Missing,Float64},
            Union{Missing,Float64},
        },
    }[]
    for curve in result.segmented_curves
        per_distribution = get(result.total_average_R_T, curve.label, nothing)
        per_distribution === nothing && continue
        for distribution in result.mass_yields
            entry = get(per_distribution, distribution.label, nothing)
            entry === nothing && continue
            more = if _is_trend(curve)
                get(result.segment_count_sensitivity, distribution.label, Tuple{Int,Float64}[])
            else
                Tuple{Int,Float64}[]
            end
            alternative(extra) = begin
                found = findfirst(t -> t[1] == segment_count(curve.fit) + extra, more)
                found === nothing ? missing : rounded(more[found][2])
            end
            spread =
                _is_trend(curve) ? leave_one_out_spread(result, distribution.label) : nothing
            spread_field(field) = spread === nothing ? missing : rounded(spread[field])
            push!(
                rows,
                (
                    segmented_curve = curve.label,
                    mass_yield = distribution.label,
                    R_T = rounded(entry.value),
                    R_T_uncertainty = rounded(entry.uncertainty),
                    R_T_uncertainty_independent_points = rounded(
                        entry.uncertainty_independent_points
                    ),
                    R_T_uncertainty_leave_one_out = spread_field(:uncertainty),
                    R_T_leave_one_out_min = spread_field(:min),
                    R_T_leave_one_out_max = spread_field(:max),
                    yield_fraction = rounded(entry.yield_fraction),
                    R_T_one_more_segment = alternative(1),
                    R_T_two_more_segments = alternative(2),
                ),
            )
        end
    end
    return rows
end

# The systematic trend refitted with each pooled dataset left out: one row per dataset and yield
# distribution it has ⟨R_T⟩ over, a single row without a distribution where the refit has none,
# so that every refit, a failed one included, is on record.
function _leave_one_out_rows(
    result::ExtractionResult, accessions::AbstractDict, digits::Integer
)
    rows = NamedTuple{
        (
            :dataset_left_out,
            :accession,
            :mass_yield,
            :R_T,
            :segments,
            :breakpoints,
            :reduced_chi_squared,
            :outcome,
        ),
        Tuple{
            String,
            Union{Missing,String},
            Union{Missing,String},
            Union{Missing,Float64},
            Union{Missing,Int},
            Union{Missing,String},
            Union{Missing,Float64},
            String,
        },
    }[]
    for entry in result.leave_one_out
        common = (
            dataset_left_out = join(entry.datasets, " + "),
            accession = _text_field(
                join(filter(!isempty, [get(accessions, l, "") for l in entry.datasets]), " ")
            ),
        )
        refit = (
            segments = entry.segments == 0 ? missing : entry.segments,
            breakpoints = _text_field(join(entry.breakpoints, " ")),
            reduced_chi_squared = if isnan(entry.reduced_chi_squared)
                missing
            else
                round(entry.reduced_chi_squared; sigdigits = digits)
            end,
            outcome = entry.outcome,
        )
        if isempty(entry.total_average_R_T)
            push!(rows, (; common..., mass_yield = missing, R_T = missing, refit...))
            continue
        end
        for distribution in result.mass_yields
            value = get(entry.total_average_R_T, distribution.label, nothing)
            value === nothing && continue
            push!(
                rows,
                (;
                    common...,
                    mass_yield = distribution.label,
                    R_T = round(value; sigdigits = digits),
                    refit...,
                ),
            )
        end
    end
    return rows
end

# Per-dataset diagnostics, written for every dataset whether or not it was pooled or fitted, with
# the reaction-code qualifiers the retrieval recorded for it.
function _diagnostic_rows(result::ExtractionResult, accessions::AbstractDict, digits::Integer)
    excluded = _excluded_labels(result.configuration.excluded_datasets, result.datasets)
    pooled = Set(pooled_datasets(result))
    superseded = superseded_datasets(result)
    # The other datasets of its experiment, for a dataset pooled as one with them.
    partners = Dict{String,Vector{String}}(
        label => filter(!=(label), group) for group in result.correlation_groups for
        label in group
    )
    return [
        (
            dataset = d.label,
            accession = _text_field(accessions[d.label]),
            points = d.points,
            pairs = d.pairs,
            first_pair = d.first_pair,
            last_pair = d.last_pair,
            coverage = round(d.coverage; sigdigits = digits),
            outside_physical_range = d.outside_physical_range,
            complement_sum = round(d.complement_sum; sigdigits = digits),
            complement_spread = round(d.complement_spread; sigdigits = digits),
            without_uncertainties = d.without_uncertainties,
            qualifiers = _text_field(join(_qualifiers(data.record), " ")),
            pooling_weight = round(_measured_fraction(data.record); sigdigits = digits),
            flagged = !isempty(_flagged(_qualifiers(data.record))),
            pair_sum_deviation = _scale_field(data, :deviation, digits),
            scale_consistent = _scale_field(data, :consistent, digits),
            pooled = d.label in pooled,
            pooled_with = _text_field(join(get(partners, d.label, String[]), " + ")),
            exclusion_reason = _text_field(
                get(
                    excluded,
                    d.label,
                    haskey(superseded, d.label) ? "superseded by $(superseded[d.label])" : "",
                ),
            ),
            segmented_curve = _text_field(get(result.dataset_outcomes, d.label, "")),
            curve_flagged = haskey(result.curve_flags, d.label),
            curve_flag_reason = _text_field(get(result.curve_flags, d.label, "")),
        ) for (d, data) in zip(result.dataset_diagnostics, result.datasets)
    ]
end

# The longest file name, in bytes, that common file systems admit (NAME_MAX).
const MAX_FILE_NAME_BYTES = 255

# The files of a run directory that carry the run identifier in their names.
function _identifier_file_names(identifier::AbstractString)
    return (
        "manifest_$(identifier).toml",
        "segmented_curves_$(identifier).csv",
        "total_average_R_T_$(identifier).csv",
        "leave_one_out_$(identifier).csv",
    )
end

# Every file name a run writes: those named by the identifier and the tables of each curve.
function _file_names(identifier::AbstractString, tokens::AbstractDict)
    names = collect(String, _identifier_file_names(identifier))
    for token in values(tokens)
        push!(
            names,
            "r_nu_vs_A_H_$(token).csv",
            "R_T_vs_A_H_$(token).csv",
            "r_nu_vs_A_H_segmented_$(token).csv",
            "R_T_vs_A_H_segmented_$(token).csv",
            "r_nu_vs_A_H_pivots_$(token).csv",
        )
    end
    return names
end

_file_token(label::AbstractString) = replace(strip(label), r"[^A-Za-z0-9.\-]+" => "_")

# A text field of a table: one that holds nothing is written as an empty field, as an absent
# number is, and not as a quoted empty string.
_text_field(text::AbstractString) = isempty(text) ? missing : String(text)
