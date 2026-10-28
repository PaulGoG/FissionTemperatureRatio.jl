# Run identification and metadata, so that every result traces back to a configuration, a commit
# and a machine.

"""
    RUN_IDENTIFIER_ABBREVIATIONS

The short key each configuration key is recorded under in [`run_parameters`](@ref), keyed by the
key's path in the configuration file.

The run identifier is formed from these short keys: the few behind the fragmentation domain are
readable tokens of it, and every one enters its hash `h` through the canonical text of the
parameters (see [`run_identifier`](@ref)). Every key [`run_parameters`](@ref) uses has an entry,
which the test suite asserts; a key without one is an error rather than a silently invented
short key.
"""
const RUN_IDENTIFIER_ABBREVIATIONS = Dict(
    "system.incident_energy" => "E",
    "fragmentation.charges_per_mass" => "nZ",
    "fragmentation.heavy_mass_min" => "AHmin",
    "fragmentation.heavy_mass_max" => "AHmax",
    "fragmentation.charge_distribution_file" => "chg",
    "fragmentation.zero_polarization_at_symmetry" => "dZ0",
    "level_density.model" => "ldm",
    "level_density.mass_excess_file" => "mass",
    "level_density.shell_correction_file" => "sc",
    "level_density.deformed_branch" => "def",
    "level_density.ratio_averaging" => "avg",
    "level_density.mean_kinetic_energy_file" => "TKE",
    "multiplicity.exclude" => "excl",
    "yield.subdirectory" => "Y",
    "yield.mass_yield_file" => "Yf",
    "yield.symmetrize" => "Ysym",
    "yield.min_yield_coverage" => "Ycov",
    "segments.max_segments" => "nseg",
    "segments.min_points_per_segment" => "npts",
    "segments.pin_symmetric_split" => "pin",
    "segments.required_windows" => "win",
    "segments.windows_apply_to_datasets" => "wdat",
    "segments.min_pair_coverage" => "cov",
    "segments.autocorrelation_lags" => "lags",
)

# Keys that change the result of a Gilbert-Cameron run only; a back-shifted Fermi gas run records
# neither, so neither enters its hash.
const GILBERT_CAMERON_KEYS = (
    "level_density.shell_correction_file", "level_density.deformed_branch"
)

# The keys that stay readable in the identifier: those behind the [domain] table of the manifest,
# which a consuming code checks a curve against, the zero polarization at the symmetric split
# among them.
const READABLE_IDENTIFIER_KEYS = (
    "level_density.model",
    "level_density.ratio_averaging",
    "fragmentation.charges_per_mass",
    "fragmentation.charge_distribution_file",
    "level_density.mass_excess_file",
    "fragmentation.zero_polarization_at_symmetry",
    "level_density.deformed_branch",
)

# Readable tokens that name an input: a named source by its name, a path by its hash.
const READABLE_INPUT_TOKENS = ("chg", "mass")

const IDENTIFIER_HASH_REASON = "h is the first eight hexadecimal digits of the SHA-1 of the \
                                canonical text, every parameter on a line of its own in the \
                                order of the short keys; a path among the readable tokens is \
                                written as the same digits of the SHA-1 of its spelling"

_hash_token(canonical::AbstractString) = bytes2hex(sha1(canonical))[1:8]

# A path relative to the data directory, which is how the configuration spelt it.
function _relative_input(path::AbstractString, configuration::Configuration)
    return relpath(path, configuration.data_directory)
end

"""
    run_parameters(configuration) -> Dict{String,Any}

The configuration values that change the result, in full, keyed by the short key of
[`RUN_IDENTIFIER_ABBREVIATIONS`](@ref); the run identifier is formed from them and the run
metadata records them under `[identifier.parameters]`.

A key with no entry in the abbreviation table throws rather than being given a short key on the
spot. A path is spelt relative to the data directory, so the same inputs staged on another
machine give the same parameters; a named source enters as its name, an absent input as
`"none"`. The exclusion list enters as its sorted accessions, the reasons changing nothing; the
required windows as `[first, last]` pairs. The Gilbert-Cameron branch and shell corrections
change a Gilbert-Cameron result only, and are parameters of such a run alone. A yield directory
enters as one table of the directory, the primary distribution that is its coverage reference
and the distributions excluded from it; a primary distribution alone enters under `Yf`. The
yield coverage floor is a parameter of a run that names a yield distribution, and of no other.
The system is not a parameter: it names the directory the identifier sits in.
`significant_digits` changes how a number is rendered, not the number, and is left out.

Every value is one TOML represents, so the table read back from `metadata.toml` gives the same
identifier through [`run_identifier`](@ref).
"""
function run_parameters(configuration::Configuration)
    fragmentation = configuration.fragmentation
    level_density = configuration.level_density
    segments = configuration.segments
    gilbert_cameron = level_density.model == "GC"
    entries = Pair{String,Any}[
        "system.incident_energy" => configuration.system.incident_energy,
        "fragmentation.charges_per_mass" => fragmentation.charges_per_mass,
        "fragmentation.heavy_mass_min" => fragmentation.heavy_mass_min,
        "fragmentation.heavy_mass_max" => fragmentation.heavy_mass_max,
        "fragmentation.charge_distribution_file" => _input_value(
            fragmentation.charge_distribution, configuration
        ),
        "fragmentation.zero_polarization_at_symmetry" => fragmentation.zero_polarization_at_symmetry,
        "level_density.model" => level_density.model,
        "level_density.mass_excess_file" => _input_value(
            level_density.mass_excess_file, configuration
        ),
        "level_density.shell_correction_file" => _input_value(
            level_density.shell_correction_file, configuration
        ),
        "level_density.deformed_branch" => level_density.deformed_branch,
        "level_density.ratio_averaging" => level_density.ratio_averaging,
        "level_density.mean_kinetic_energy_file" => _input_value(
            level_density.mean_kinetic_energy_file, configuration
        ),
        "multiplicity.exclude" => sort!(collect(keys(configuration.excluded_datasets))),
        "yield.subdirectory" => _yield_source_value(configuration),
        "yield.mass_yield_file" => _input_value(configuration.yield_file, configuration),
        "yield.symmetrize" => configuration.symmetrize_yields,
        "yield.min_yield_coverage" => configuration.min_yield_coverage,
        "segments.max_segments" => segments.max_segments,
        "segments.min_points_per_segment" => segments.min_points_per_segment,
        "segments.pin_symmetric_split" => segments.pin_symmetric_split,
        "segments.required_windows" => [[first(w), last(w)] for w in segments.required_windows],
        "segments.windows_apply_to_datasets" => segments.windows_apply_to_datasets,
        "segments.min_pair_coverage" => segments.min_pair_coverage,
        "segments.autocorrelation_lags" => segments.autocorrelation_lags,
    ]
    parameters = Dict{String,Any}()
    for (key, value) in entries
        !gilbert_cameron && key in GILBERT_CAMERON_KEYS && continue
        # The primary distribution alone enters as `Yf`; a directory enters as `Y`, one table of
        # the directory, the primary that is its coverage reference and the exclusions.
        directory = configuration.yield_directory !== nothing
        key == "yield.subdirectory" &&
            !directory &&
            configuration.yield_file !== nothing &&
            continue
        key == "yield.mass_yield_file" &&
            (directory || configuration.yield_file === nothing) &&
            continue
        key == "yield.min_yield_coverage" && configuration.yield_file === nothing && continue
        haskey(RUN_IDENTIFIER_ABBREVIATIONS, key) || throw(
            ArgumentError("no entry in RUN_IDENTIFIER_ABBREVIATIONS for the key $(repr(key))"),
        )
        parameters[RUN_IDENTIFIER_ABBREVIATIONS[key]] = value
    end
    return parameters
end

# A configured input: a named source ("wahl", "mean", "ame2020", "gilbert_cameron_1965") by its
# name, a path by its spelling relative to the data directory, an absent one as "none".
const NAMED_SOURCES = ("wahl", "mean", SHIPPED_MASS_TABLE, SHIPPED_SHELL_CORRECTIONS)

_input_value(::Nothing, ::Configuration) = "none"
function _input_value(value::AbstractString, configuration::Configuration)
    value in NAMED_SOURCES && return String(value)
    return _relative_input(value, configuration)
end

# A directory of yield distributions, with the primary that is its coverage reference and the
# distributions excluded from it, as one parameter; "none" where no directory is read.
function _yield_source_value(configuration::Configuration)
    configuration.yield_directory === nothing &&
        return _input_value(configuration.yield_directory, configuration)
    return Dict{String,Any}(
        "subdirectory" => _input_value(configuration.yield_directory, configuration),
        "mass_yield_file" => _input_value(configuration.yield_file, configuration),
        "exclude" => sort!(collect(keys(configuration.excluded_mass_yields))),
    )
end

# The yield distributions a run read but took no total average over, each with the reason.
function _not_averaged(result::ExtractionResult)
    configuration = result.configuration
    floor = configuration.min_yield_coverage
    excluded = if configuration.yield_directory === nothing
        Dict{String,String}()
    else
        _excluded_labels(
            configuration.excluded_mass_yields,
            result.mass_yields,
            result.mass_yield_records,
        )
    end
    reasons = Dict{String,Any}()
    for (label, coverage) in result.mass_yield_coverage
        if haskey(excluded, label)
            reasons[label] = "excluded by configuration: $(excluded[label])"
        elseif coverage < floor
            reasons[label] = "coverage $(round(coverage; digits = 3)) below the floor $(floor)"
        end
    end
    return reasons
end

"""
    run_identifier(configuration) -> String
    run_identifier(parameters::AbstractDict) -> String

The token that names one run of one fissioning system: DrWatson's `savename` over the readable
tokens and `h`, sorted and joined by `_`.

The readable tokens are the parameters behind the `[domain]` table of the manifest, which a
consuming code checks a curve against — `ldm`, `avg`, `nZ`, `chg`, `mass`, `dZ0`, and `def` for
a Gilbert-Cameron run. `dZ0` is the configured zero polarization at the symmetric split; the
table records the one the domain applied, which is the same, since the configuration refuses the
choice with Wahl's model, the one charge model that would ignore it. A named charge distribution
or mass table appears by its name, a path as
the first eight hexadecimal digits of the SHA-1 of its spelling. `h` is the first eight
hexadecimal digits of the SHA-1 of the canonical text of the whole of
[`run_parameters`](@ref), readable tokens included: one `key=value` line per parameter in the
order of the keys, strings quoted, integer-valued reals without a decimal part, lists in
brackets and tables in braces. Every other parameter therefore enters the identifier through
`h` alone.

Two configurations that differ in any result-changing key give different identifiers, and two
that agree in all of them give the same one. The multiplicity directory is not a parameter: it
is `<system>/nu_vs_A` by convention, so two runs of one system on different directories share
an identifier. The second method takes the parameters as
`metadata.toml` records them under `[identifier.parameters]`, so the identifier of a run is
recomputed from its own record. The run directory is `data/sims/<system>/<identifier>/`.

# Examples
```julia
run_identifier(configuration)
# "avg=charge_resolved_chg=wahl_dZ0=false_h=…_ldm=BSFG_mass=ame2020_nZ=5"
```
"""
run_identifier(configuration::Configuration) = run_identifier(run_parameters(configuration))

function run_identifier(parameters::AbstractDict)
    return savename(_identifier_tokens(parameters); connector = "_", val_to_string = _token)
end

# The tokens of the run identifier: the readable parameters, a path among them as its hash, and
# h. A missing readable parameter is refused, `def` excepted outside a Gilbert-Cameron run.
function _identifier_tokens(parameters::AbstractDict)
    gilbert_cameron =
        get(parameters, RUN_IDENTIFIER_ABBREVIATIONS["level_density.model"], "") == "GC"
    tokens = Dict{String,Any}()
    for key in READABLE_IDENTIFIER_KEYS
        token = RUN_IDENTIFIER_ABBREVIATIONS[key]
        !gilbert_cameron && key in GILBERT_CAMERON_KEYS && continue
        haskey(parameters, token) || throw(
            ArgumentError("the parameters carry no $(repr(token)) ($(key)), a readable token \
                           of the run identifier"),
        )
        value = parameters[token]
        tokens[token] = if token in READABLE_INPUT_TOKENS && !(value in NAMED_SOURCES)
            _hash_token(value)
        else
            value
        end
    end
    tokens["h"] = _hash_token(_canonical_parameters(parameters))
    return tokens
end

# The canonical text h is taken over. It depends on the values alone, not on how they were read,
# so the parameters of a configuration and the same parameters read back from metadata.toml give
# the same text: strings quoted, integer-valued reals without a decimal part, other reals in their
# shortest round-trip form.
function _canonical_parameters(parameters::AbstractDict)
    return join(
        (
            "$(key)=$(_canonical_value(parameters[key]))" for
            key in sort!(collect(keys(parameters)))
        ),
        "\n",
    )
end

_canonical_value(x::Bool) = string(x)
_canonical_value(x::Real) = _token(x)
_canonical_value(x::AbstractString) = repr(String(x))
_canonical_value(x::AbstractVector) = "[" * join((_canonical_value(v) for v in x), ",") * "]"
function _canonical_value(x::AbstractDict)
    entries = ("$(key)=$(_canonical_value(x[key]))" for key in sort!(collect(keys(x))))
    return "{" * join(entries, ",") * "}"
end

_token(x::Bool) = string(x)
_token(x::Integer) = string(x)
_token(x::AbstractString) = String(x)
# Trailing zeros carry no information and would make two spellings of one value.
_token(x::Real) = isinteger(x) ? string(Int(round(x))) : string(x)

"""
    run_metadata(result) -> Dict{String,Any}

What the library knows about a run, as nested tables ready for `TOML.print`: the system, every
configuration value in full, the identifier tokens with the parameters and the canonical text
they were formed from, the
fragmentation domain as the manifest records it together with the charge model that built it,
the input files and the retrievals that produced them — the parser revision and the record of the
`⟨TKE⟩(A)` dataset in particular, with the heavy masses it was interpolated or extrapolated to —
the package and dependency versions, and a summary of the result: each segmented curve, the
identities checked at the symmetric split, what became of every dataset, the autocorrelation the
covariance of the trend was built with, the leave-one-out refits of the trend, and the flagged
dataset curves.

The caller adds what only it knows — the run identifier's place on disk, the configuration path,
the files written, the commit, the platform — and writes the whole; `scripts/run.jl` does.
"""
function run_metadata(result::ExtractionResult)
    configuration = result.configuration
    fragmentation = configuration.fragmentation
    level_density = configuration.level_density
    segments = configuration.segments
    relative(path) = path === nothing ? "none" : _relative_input(path, configuration)
    parameters = run_parameters(configuration)

    curves = Dict{String,Any}(
        curve.label => Dict{String,Any}(
            "kind" => _kind(curve),
            "segments" => FissionTemperatureRatio.segments(curve.fit),
            "breakpoints" => curve.fit.breakpoints,
            "pivots" => [[point[1], point[2]] for point in pivots(curve.fit)],
            "reduced_chi_squared" => curve.fit.wrss / curve.fit.dof,
            "chi_squared_over_expectation" => curve.fit.wrss / curve.fit.expected_wrss,
            "correlated_points" => curve.fit.correlated,
            "bic" => curve.fit.bic,
            "bic_by_order" => [[order, value] for (order, value) in curve.fit.selection],
            "weights_imputed" => curve.fit.weights_imputed,
            "points" => curve.fit.points,
            "measured_points" => curve.fit.measured_points,
            "pairs" => curve.fit.points,
            "coverage" => curve.coverage,
            "range_mean_R_T" => collect(result.range_mean_R_T[curve.label]),
        ) for curve in result.segmented_curves
    )
    domain = manifest_domain(result)
    symmetry = result.symmetry
    accessions = curve_accessions(result)
    flagged = Dict{String,Any}(
        data.label => _qualifiers(data.record) for
        data in result.datasets if !isempty(_flagged(_qualifiers(data.record)))
    )
    return Dict{String,Any}(
        "system" => system_record(configuration.system),
        "configuration" => Dict{String,Any}(
            "charges_per_mass" => fragmentation.charges_per_mass,
            "heavy_mass_min" => fragmentation.heavy_mass_min,
            "heavy_mass_max" => fragmentation.heavy_mass_max,
            "charge_distribution_file" =>
                _input_value(fragmentation.charge_distribution, configuration),
            "zero_polarization_at_symmetry" => fragmentation.zero_polarization_at_symmetry,
            "model" => level_density.model,
            "mass_excess_file" =>
                _input_value(level_density.mass_excess_file, configuration),
            "deformed_branch" => level_density.deformed_branch,
            "ratio_averaging" => level_density.ratio_averaging,
            "mean_kinetic_energy_file" =>
                _input_value(level_density.mean_kinetic_energy_file, configuration),
            "excluded_datasets" => Dict{String,Any}(configuration.excluded_datasets),
            "excluded_mass_yields" => Dict{String,Any}(configuration.excluded_mass_yields),
            "max_segments" => segments.max_segments,
            "min_points_per_segment" => segments.min_points_per_segment,
            "pin_symmetric_split" => segments.pin_symmetric_split,
            "required_windows" => [[first(w), last(w)] for w in segments.required_windows],
            "windows_apply_to_datasets" => segments.windows_apply_to_datasets,
            "min_pair_coverage" => segments.min_pair_coverage,
            "autocorrelation_lags" => segments.autocorrelation_lags,
            "symmetrize_yields" => configuration.symmetrize_yields,
            "min_yield_coverage" => configuration.min_yield_coverage,
            "min_package_version" => string(configuration.min_retrieval_version),
            "significant_digits" => configuration.output.significant_digits,
        ),
        # The identifier's tokens, every parameter it was formed from, and the text h was taken
        # over, so that the identifier is recomputed from this record alone.
        "identifier" => Dict{String,Any}(
            "tokens" => _identifier_tokens(parameters),
            "parameters" => parameters,
            "canonical" => _canonical_parameters(parameters),
            "hash_reason" => IDENTIFIER_HASH_REASON,
        ),
        # The [domain] table of the manifest, field for field, and what it abbreviates.
        "domain" => Dict{String,Any}(
            "level_density_model" => domain.level_density_model,
            "deformed_branch" => domain.deformed_branch,
            "ratio_averaging" => domain.ratio_averaging,
            "excitation_weighted" => domain.excitation_weighted,
            "charges_per_mass" => domain.charges_per_mass,
            "charge_model" => domain.charge_model,
            "mass_table" => domain.mass_table,
            "package_version" => domain.package_version,
            "zero_polarization_at_symmetry" => domain.zero_polarization_at_symmetry,
            "fragmentations" => length(result.domain),
        ),
        "inputs" => Dict{String,Any}(
            "mass_excess_file" =>
                _input_value(level_density.mass_excess_file, configuration),
            "shell_correction_file" =>
                _input_value(level_density.shell_correction_file, configuration),
            "charge_distribution_file" =>
                _input_value(fragmentation.charge_distribution, configuration),
            "multiplicity_directory" => relative(configuration.multiplicity_directory),
            "yield_directory" => relative(configuration.yield_directory),
            "mass_yield_file" => relative(configuration.yield_file),
            # The EXFOR accession of every yield distribution averaged over, from its retrieval
            # record.
            "mass_yield_accessions" => Dict{String,Any}(
                y.label => _accession(
                    get(result.mass_yield_records, y.label, nothing), y.source
                ) for y in result.mass_yields
            ),
            "datasets" => [basename(data.source) for data in result.datasets],
            "mass_yields" => [basename(data.source) for data in result.mass_yields],
            "flagged_datasets" => flagged,
            # Interpolated datasets, pooled at raw points over rows written.
            "pooling_weights" => Dict{String,Any}(
                data.label => w for data in result.datasets for
                w in (_measured_fraction(data.record),) if w < 1
            ),
            # Datasets whose pair sum the retrieval finds off the scale of ν̄: used, not
            # corrected, since a uniform scale cancels in r_ν. [deviation, uncertainty].
            "scale_inconsistent_datasets" => Dict{String,Any}(
                data.label => [s.deviation, coalesce(s.uncertainty, NaN)] for
                data in result.datasets for s in (_pair_sum(data.record),) if
                s !== nothing && s.consistent === false
            ),
            "flagged_mass_yields" => Dict{String,Any}(
                label => record.qualifiers for
                (label, record) in result.mass_yield_records if
                !isempty(_flagged(record.qualifiers))
            ),
            "retrievals" => _retrieval_runs(configuration),
        ),
        "mean_kinetic_energy" => _kinetic_energy_record(result),
        "source" => Dict{String,Any}(
            "package_version" => string(pkgversion(@__MODULE__)),
            "dependencies" => _dependency_versions(),
        ),
        "result" => Dict{String,Any}(
            "segmented_curves" => curves,
            "dataset_outcomes" => Dict{String,Any}(result.dataset_outcomes),
            # The share of the primary distribution's heavy-fragment yield each distribution
            # covers, and the ones not averaged over, with the reason.
            "mass_yield_coverage" => Dict{String,Any}(result.mass_yield_coverage),
            "mass_yields_not_averaged" => _not_averaged(result),
            # How the covariance of the systematic trend was built: the correlogram of the pooled
            # datasets' deviations, the autocorrelation fitted to it, the scale of the covariance
            # against the expectation of χ², and the leave-one-out spread of ⟨R_T⟩ per yield.
            "trend_uncertainty" => _trend_uncertainty_record(result),
            # The pooled datasets that belong to one experiment and were combined into one before
            # pooling, with the relation their retrieval record states.
            "correlation_groups" => [
                _correlation_group_record(result, group) for group in result.correlation_groups
            ],
            # Datasets a republication supersedes, each with its successor, through which alone
            # it enters the pool.
            "superseded_datasets" => Dict{String,Any}(superseded_datasets(result)),
            # Dataset curves that resolve no minimum, with the reason; they stay in the manifest.
            "flagged_curves" => Dict{String,Any}(result.curve_flags),
            # Whether the refits with each pooled dataset left out were asked for, and the trend
            # so refitted, in the order read.
            "leave_one_out_requested" => result.leave_one_out_requested,
            "leave_one_out" => [
                _leave_one_out_record(entry, accessions) for entry in result.leave_one_out
            ],
            # ⟨R_T⟩ of the systematic trend at the selected number of segments and one and two
            # more, per yield distribution: [segments, ⟨R_T⟩].
            "segment_count_sensitivity" => Dict{String,Any}(
                label => [[k, value] for (k, value) in entries] for
                (label, entries) in result.segment_count_sensitivity
            ),
            "symmetry" => Dict{String,Any}(
                "charge_set_invariant" => if symmetry.charge_set_invariant === missing
                    "not applicable"
                else
                    symmetry.charge_set_invariant
                end,
                "R_T_at_symmetric_split" => if symmetry.R_T_at_symmetric_split === missing
                        "not applicable"
                    else
                        symmetry.R_T_at_symmetric_split
                    end,
                "pinned_curves" => symmetry.pinned_curves,
                "warnings" => symmetry.warnings,
            ),
        ),
    )
end

# The correlogram, with NaN at a lag that holds no value, whether an autocorrelation could be
# fitted to it and the coefficient the covariance was formed with, zero where none could, the
# scale of the trend's covariance against the expectation of χ², and the leave-one-out spread of
# ⟨R_T⟩ over each yield distribution it was formed for.
function _trend_uncertainty_record(result::ExtractionResult)
    ρ = deviation_autocorrelation(result)
    record = Dict{String,Any}(
        "deviation_correlogram" =>
            Float64[coalesce(c, NaN) for c in result.deviation_correlogram],
        "autocorrelation_lags" => result.configuration.segments.autocorrelation_lags,
        "autocorrelation_estimated" => !ismissing(ρ),
        "deviation_autocorrelation" => coalesce(ρ, 0.0),
    )
    index = findfirst(c -> c.label == SYSTEMATIC_TREND_LABEL, result.segmented_curves)
    if index !== nothing
        fit = result.segmented_curves[index].fit
        record["reduced_chi_squared"] = fit.wrss / fit.dof
        record["expected_wrss"] = fit.expected_wrss
        record["chi_squared_over_expectation"] = fit.wrss / fit.expected_wrss
        record["covariance_scale"] = _covariance_scale(
            fit.wrss, fit.expected_wrss, fit.weights_imputed, fit.points
        )
    end
    spreads = Dict{String,Any}()
    for distribution in result.mass_yields
        spread = leave_one_out_spread(result, distribution.label)
        spread === nothing && continue
        spreads[distribution.label] = Dict{String,Any}(
            "min" => spread.min,
            "max" => spread.max,
            "uncertainty" => spread.uncertainty,
            "refits" => spread.refits,
        )
    end
    record["leave_one_out"] = spreads
    return record
end

# One experiment pooled as one: its datasets, their accessions, the relation every record
# states, which decides how its curve was formed, and what each record states.
function _correlation_group_record(result::ExtractionResult, group::Vector{String})
    position = Dict(data.label => index for (index, data) in enumerate(result.datasets))
    indices = Int[position[label] for label in group]
    relation = _relation(result.datasets, indices)
    return Dict{String,Any}(
        "datasets" => group,
        "accessions" => [_accession(result.datasets[i]) for i in indices],
        "relation" => isempty(relation) ? "unstated or differing" : relation,
        "relations_stated" =>
            [_correlation_relation(result.datasets[i].record) for i in indices],
    )
end

# One leave-one-out refit of the trend; the reduced chi-squared is omitted where there is no fit.
function _leave_one_out_record(entry::LeaveOneOut, accessions::AbstractDict)
    record = Dict{String,Any}(
        "datasets" => entry.datasets,
        "accessions" => [get(accessions, label, "") for label in entry.datasets],
        "outcome" => entry.outcome,
        "segments" => entry.segments,
        "breakpoints" => entry.breakpoints,
        "total_average_R_T" => Dict{String,Any}(entry.total_average_R_T),
    )
    isnan(entry.reduced_chi_squared) ||
        (record["reduced_chi_squared"] = entry.reduced_chi_squared)
    return record
end

# The `[run]` table of every retrieval record in the input directories: which parser revision
# produced the files, from which of its configurations, and when.
function _retrieval_runs(configuration::Configuration)
    runs = Dict{String,Any}()
    for directory in (
        configuration.multiplicity_directory,
        configuration.yield_directory,
        _parent(configuration.yield_file),
        _parent(configuration.level_density.mean_kinetic_energy_file),
    )
        directory === nothing && continue
        for (record, document) in _retrieval_documents(directory)
            run = get(document, "run", nothing)
            runs[_relative_input(record, configuration)] =
                run isa AbstractDict ? Dict{String,Any}(run) : "no [run] table"
        end
    end
    return runs
end

function _offset_record(energies, distribution, system)
    offset = mean_kinetic_energy_offset(energies, distribution, system)
    # The heavy masses the distribution gives: a partial one biases the mean towards its own.
    record = Dict{String,Any}(
        "heavy_masses" => length(offset.masses),
        "first_A_H" => first(offset.masses),
        "last_A_H" => last(offset.masses),
        "mean_MeV" => offset.mean,
    )
    offset.standard === nothing &&
        return merge(record, Dict{String,Any}("standard" => "none recorded for this system"))
    return merge(
        record,
        Dict{String,Any}(
            "standard_MeV" => offset.standard,
            "standard_uncertainty_MeV" => offset.standard_uncertainty,
            "offset_MeV" => offset.offset,
        ),
    )
end

_parent(::Nothing) = nothing
_parent(path::AbstractString) = dirname(path)

# The ⟨TKE⟩(A) dataset that weighted the inversion: the file, its retrieval record entry, the
# parser revision that wrote it, a digest of the whole record, and the heavy masses that rest on
# a neighbour's value rather than their own.
function _kinetic_energy_record(result::ExtractionResult)
    configuration = result.configuration
    energies = result.mean_kinetic_energy
    weighted = manifest_domain(result).excitation_weighted
    energies === nothing && return Dict{String,Any}(
        "excitation_weighted" => weighted,
        "file" => "none",
        "reason" => if result.averaging isa ChargeResolved
            "no ⟨TKE⟩(A) dataset configured; weights p(Z, A_H) alone"
        else
            "ratio_averaging = $(repr(configuration.level_density.ratio_averaging)) forms an \
             effective ratio and takes no excitation weight"
        end,
    )
    entry = Dict{String,Any}(
        "excitation_weighted" => weighted,
        "file" => _relative_input(energies.source, configuration),
        "accession" => _accession(energies),
        "label" => energies.label,
        "unit" => "MeV",
        "measured" => energies.measured,
        "averaged_with_complement" => energies.averaged,
        "interpolated" => energies.interpolated,
        "extrapolated" => energies.extrapolated,
        # The input's yield-weighted mean against the energy standard, over each distribution
        # the run read: the scale a double-energy measurement sits on.
        "offset_from_standard" => Dict{String,Any}(
            distribution.label =>
                _offset_record(energies, distribution, configuration.system) for
            distribution in result.mass_yields
        ),
    )
    record = energies.record
    if record !== nothing
        entry["retrieval_record"] = _relative_input(record.path, configuration)
        entry["retrieval_record_sha1"] = bytes2hex(open(sha1, record.path))
        entry["retrieval_run"] = record.run
        entry["retrieval_entry"] = record.entry
        entry["qualifiers"] = record.qualifiers
        entry["flagged"] = _flagged(record.qualifiers)
    end
    return entry
end

# The versions the run used, read from the active manifest; manifests are not version-controlled.
function _dependency_versions()
    project = Base.active_project()
    project === nothing && return Dict{String,Any}("status" => "unavailable")
    manifest = joinpath(dirname(project), "Manifest.toml")
    (isfile(project) && isfile(manifest)) || return Dict{String,Any}("status" => "unavailable")

    versions = Dict{String,Any}()
    try
        direct = keys(get(TOML.parsefile(project), "deps", Dict{String,Any}()))
        document = TOML.parsefile(manifest)
        # A manifest of format 2 nests its entries under `deps`; earlier ones list them at the top
        # level alongside a few scalar keys.
        entries = get(document, "deps") do
            return Dict(k => v for (k, v) in document if v isa AbstractVector)
        end
        for name in direct
            haskey(entries, name) || continue
            for entry in entries[name]
                version = get(entry, "version", nothing)
                version === nothing || (versions[name] = version)
            end
        end
    catch exception
        exception isa TOML.ParserError || rethrow()
        return Dict{String,Any}("status" => "unavailable")
    end
    return versions
end

# Append a numeric suffix rather than overwrite, in the manner of DrWatson's safesave, so that a
# figure written twice never destroys the first.
function _unused_path(path::AbstractString)
    isfile(path) || return String(path)
    stem, extension = splitext(path)
    index = 1
    while isfile("$(stem)#$(index)$(extension)")
        index += 1
    end
    return "$(stem)#$(index)$(extension)"
end
