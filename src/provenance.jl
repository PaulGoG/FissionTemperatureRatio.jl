# Run identification and metadata, so that every result traces back to a configuration, a commit
# and a machine.

# A path relative to the data directory, which is how the configuration spelt it.
function _relative_input(path::AbstractString, configuration::Configuration)
    return relpath(path, configuration.data_directory)
end

# A configured input: a named source ("wahl", "mean", "ame2020", "gilbert_cameron_1965") by its
# name, a path by its spelling relative to the data directory.
const NAMED_SOURCES = ("wahl", "mean", SHIPPED_MASS_TABLE, SHIPPED_SHELL_CORRECTIONS)

function _input_value(value::AbstractString, configuration::Configuration)
    value in NAMED_SOURCES && return String(value)
    return _relative_input(value, configuration)
end

# Every value a configuration holds, defaults included, by section and key as the file spells
# them; an optional input that is not given, and the [yield] section of a run without one, are
# absent. Paths are relative to the data directory.
function _configuration_values(configuration::Configuration)
    fragmentation = configuration.fragmentation
    level_density = configuration.level_density
    segments = configuration.segments
    input(value) = _input_value(value, configuration)
    sections = Dict{String,Any}(
        "system" => Dict{String,Any}(
            "target_A" => configuration.system.target.A,
            "target_Z" => configuration.system.target.Z,
            "channel" => configuration.system.channel,
            "incident_energy" => configuration.system.incident_energy,
        ),
        "fragmentation" => Dict{String,Any}(
            "charges_per_mass" => fragmentation.charges_per_mass,
            "heavy_mass_min" => fragmentation.heavy_mass_min,
            "heavy_mass_max" => fragmentation.heavy_mass_max,
            "charge_distribution_file" => input(fragmentation.charge_distribution),
            "zero_polarization_at_symmetry" => fragmentation.zero_polarization_at_symmetry,
        ),
        "level_density" => Dict{String,Any}(
            "model" => level_density.model,
            "mass_excess_file" => input(level_density.mass_excess_file),
            "shell_correction_file" => input(level_density.shell_correction_file),
            "deformed_branch" => level_density.deformed_branch,
            "ratio_averaging" => level_density.ratio_averaging,
        ),
        "multiplicity" => Dict{String,Any}(
            "subdirectory" =>
                _relative_input(configuration.multiplicity_directory, configuration),
            "exclude" => Dict{String,Any}(configuration.excluded_datasets),
        ),
        "segments" => Dict{String,Any}(
            "max_segments" => segments.max_segments,
            "min_points_per_segment" => segments.min_points_per_segment,
            "pin_symmetric_split" => segments.pin_symmetric_split,
            "required_windows" => [[first(w), last(w)] for w in segments.required_windows],
            "windows_apply_to_datasets" => segments.windows_apply_to_datasets,
            "min_pair_coverage" => segments.min_pair_coverage,
            "autocorrelation_lags" => segments.autocorrelation_lags,
        ),
        "retrieval" => Dict{String,Any}(
            "min_package_version" => string(configuration.min_retrieval_version)
        ),
        "output" =>
            Dict{String,Any}("significant_digits" => configuration.output.significant_digits),
    )
    energies = level_density.mean_kinetic_energy_file
    energies === nothing ||
        (sections["level_density"]["mean_kinetic_energy_file"] = input(energies))
    if configuration.yield_file !== nothing
        yields = Dict{String,Any}(
            "mass_yield_file" => input(configuration.yield_file),
            "symmetrize" => configuration.symmetrize_yields,
            "min_yield_coverage" => configuration.min_yield_coverage,
        )
        if configuration.yield_directory !== nothing
            yields["subdirectory"] = input(configuration.yield_directory)
            yields["exclude"] = Dict{String,Any}(configuration.excluded_mass_yields)
        end
        sections["yield"] = yields
    end
    return sections
end

"""
    run_parameters(configuration) -> Dict{String,Any}

The configuration values that change the result, in full, by section and key as the
configuration file spells them; the run identifier is formed from them and the run metadata
records them under `[identifier.parameters]`.

A value enters as the configuration holds it once read, defaults included, and a path relative
to the data directory, so the same inputs staged on another machine give the same parameters. A
key whose value cannot change this run's result is no parameter of it: the Gilbert-Cameron
branch and shell corrections of a back-shifted Fermi gas run, `windows_apply_to_datasets`
without windows, an optional input that is not given, the [yield] keys of a run without yields.
An exclusion list enters as its sorted accessions, the reasons changing nothing, and an empty
one not at all; the required windows as `[first, last]` pairs. The system is not a parameter: it
names the directory the identifier sits in. Neither are `[retrieval]`, a guard on the inputs,
nor `[output]`, which changes how a number is rendered and not the number.

Every value is one TOML represents, so the table read back from `metadata.toml` gives the same
identifier through [`run_identifier`](@ref).
"""
function run_parameters(configuration::Configuration)
    parameters = _configuration_values(configuration)
    for section in ("system", "retrieval", "output")
        delete!(parameters, section)
    end
    parameters["system"] = Dict{String,Any}(
        "incident_energy" => configuration.system.incident_energy
    )
    level_density = parameters["level_density"]
    if configuration.level_density.model != "GC"
        delete!(level_density, "shell_correction_file")
        delete!(level_density, "deformed_branch")
    end
    segments = parameters["segments"]
    if isempty(configuration.segments.required_windows)
        delete!(segments, "required_windows")
        delete!(segments, "windows_apply_to_datasets")
    end
    for section in ("multiplicity", "yield")
        haskey(parameters, section) && haskey(parameters[section], "exclude") || continue
        accessions = sort!(collect(keys(parameters[section]["exclude"])))
        if isempty(accessions)
            delete!(parameters[section], "exclude")
        else
            parameters[section]["exclude"] = accessions
        end
    end
    return parameters
end

# The readable tokens of the identifier and the parameter each stands for: those behind the
# [domain] table of the manifest, which a consuming code checks a curve against, the zero
# polarization at the symmetric split among them, and the Gilbert-Cameron branch of such a run.
const IDENTIFIER_TOKENS = (
    "avg" => ("level_density", "ratio_averaging"),
    "chg" => ("fragmentation", "charge_distribution_file"),
    "dZ0" => ("fragmentation", "zero_polarization_at_symmetry"),
    "def" => ("level_density", "deformed_branch"),
    "ldm" => ("level_density", "model"),
    "mass" => ("level_density", "mass_excess_file"),
    "nZ" => ("fragmentation", "charges_per_mass"),
)

# Tokens that name an input: a named source by its name, a path by its hash.
const IDENTIFIER_INPUT_TOKENS = ("chg", "mass")

_hash_token(canonical::AbstractString) = bytes2hex(sha1(canonical))[1:8]

"""
    run_identifier(configuration) -> String
    run_identifier(parameters::AbstractDict) -> String

The token that names one run of one fissioning system: `key=value` pairs, sorted by key and
joined by `_`, of the readable tokens and `h`.

The readable tokens are the parameters behind the `[domain]` table of the manifest, which a
consuming code checks a curve against — `ldm`, `avg`, `nZ`, `chg`, `mass`, `dZ0`, and `def` for
a Gilbert-Cameron run. `dZ0` is the configured zero polarization at the symmetric split; the
table records the one the domain applied, which is the same, since the configuration refuses the
choice with Wahl's model, the one charge model that would ignore it. A named charge distribution
or mass table appears by its name, a path as the first eight hexadecimal digits of the SHA-1 of
its spelling. `h` is the first eight hexadecimal digits of the SHA-1 of the canonical text of the
whole of [`run_parameters`](@ref), readable tokens included: one `section.key=value` line per
parameter, sorted, strings quoted, integer-valued reals without a decimal part, lists in
brackets. Every other parameter therefore enters the identifier through `h` alone.

Two configurations that differ in any result-changing key give different identifiers, and two
that agree in all of them give the same one. The second method takes the parameters as
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
    tokens = _identifier_tokens(parameters)
    return join(("$(key)=$(_token(tokens[key]))" for key in sort!(collect(keys(tokens)))), "_")
end

# The tokens of the run identifier: the readable parameters, a path among them as its hash, and
# h. A missing readable parameter is refused, `def` excepted outside a Gilbert-Cameron run.
function _identifier_tokens(parameters::AbstractDict)
    gilbert_cameron = get(get(parameters, "level_density", Dict()), "model", "") == "GC"
    tokens = Dict{String,Any}()
    for (token, (section, key)) in IDENTIFIER_TOKENS
        token == "def" && !gilbert_cameron && continue
        value = get(get(parameters, section, Dict()), key, nothing)
        value === nothing && throw(
            ArgumentError("the parameters carry no $(section).$(key), the readable token \
                           $(repr(token)) of the run identifier")
        )
        tokens[token] = if token in IDENTIFIER_INPUT_TOKENS && !(value in NAMED_SOURCES)
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
    lines = String[]
    for section in sort!(collect(keys(parameters)))
        entries = parameters[section]
        for key in sort!(collect(keys(entries)))
            push!(lines, "$(section).$(key)=$(_canonical_value(entries[key]))")
        end
    end
    return join(lines, "\n")
end

_canonical_value(x::Bool) = string(x)
_canonical_value(x::Real) = _token(x)
_canonical_value(x::AbstractString) = repr(String(x))
_canonical_value(x::AbstractVector) = "[" * join((_canonical_value(v) for v in x), ",") * "]"

_token(x::Bool) = string(x)
_token(x::Integer) = string(x)
_token(x::AbstractString) = String(x)
# Trailing zeros carry no information and would make two spellings of one value.
_token(x::Real) = isinteger(x) ? string(Int(round(x))) : string(x)

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
    run_metadata(result) -> Dict{String,Any}

What the library knows about a run, as nested tables ready for `TOML.print`: the system, every
configuration value by section and key, defaults included, the parameters the run identifier was
formed from, the input files and the retrievals that produced them — the parser revision and the
record of the `⟨TKE⟩(A)` dataset in particular, with the heavy masses it was interpolated or
extrapolated to — the package version, and a summary of the result: each segmented curve, the
identities checked at the symmetric split, what became of every dataset, the autocorrelation the
covariance of the trend was built with, the leave-one-out refits of the trend, and the flagged
dataset curves.

A key that does not apply to the run is absent rather than written with a placeholder: an
optional input not given, a test with no symmetric split to apply to, an autocorrelation no lag
could estimate. The domain the curves were extracted on is the manifest's `[domain]` table, and
the resolved dependencies are the `Manifest.toml` `scripts/run.jl` copies beside the record.

The caller adds what only it knows — the run identifier's place on disk, the configuration path,
the files written, the commit, the platform — and writes the whole; `scripts/run.jl` does.
"""
function run_metadata(result::ExtractionResult)
    configuration = result.configuration
    curves = Dict{String,Any}(
        curve.label => Dict{String,Any}(
            "kind" => _kind(curve),
            "segments" => segment_count(curve.fit),
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
            "coverage" => curve.coverage,
            "range_mean_R_T" => collect(result.range_mean_R_T[curve.label]),
        ) for curve in result.segmented_curves
    )
    accessions = curve_accessions(result)
    return Dict{String,Any}(
        "system" => system_record(configuration.system),
        "configuration" => _configuration_values(configuration),
        # Every parameter the identifier was formed from, so that it is recomputed from this
        # record alone, `run_identifier(metadata["identifier"]["parameters"])`.
        "identifier" => Dict{String,Any}("parameters" => run_parameters(configuration)),
        "inputs" => _inputs_record(result),
        "mean_kinetic_energy" => _kinetic_energy_record(result),
        "source" => Dict{String,Any}("package_version" => string(pkgversion(@__MODULE__))),
        "result" => Dict{String,Any}(
            "segmented_curves" => curves,
            "fragmentations" => length(result.domain),
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
            "symmetry" => _symmetry_record(result.symmetry),
        ),
    )
end

# The input files of a run and the retrievals that produced them.
function _inputs_record(result::ExtractionResult)
    configuration = result.configuration
    level_density = configuration.level_density
    record = Dict{String,Any}(
        "mass_excess_file" => _input_value(level_density.mass_excess_file, configuration),
        "charge_distribution_file" =>
            _input_value(configuration.fragmentation.charge_distribution, configuration),
        "multiplicity_directory" =>
            _relative_input(configuration.multiplicity_directory, configuration),
        "datasets" => [basename(data.source) for data in result.datasets],
        "mass_yields" => [basename(data.source) for data in result.mass_yields],
        # The EXFOR accession of every yield distribution read, from its retrieval record.
        "mass_yield_accessions" => Dict{String,Any}(
            y.label =>
                _accession(get(result.mass_yield_records, y.label, nothing), y.source) for
            y in result.mass_yields
        ),
        "flagged_datasets" => Dict{String,Any}(
            data.label => _qualifiers(data.record) for
            data in result.datasets if !isempty(_flagged(_qualifiers(data.record)))
        ),
        # Interpolated datasets, pooled at their measured values over their rows written.
        "pooling_weights" => Dict{String,Any}(
            data.label => w for data in result.datasets for
            w in (_measured_fraction(data.record),) if w < 1
        ),
        # Datasets whose pair sum the retrieval finds off the scale of ν̄: used, not corrected,
        # since a uniform scale cancels in r_ν. [deviation, uncertainty]; a verdict on the scale
        # is formed only where an uncertainty is stated.
        "scale_inconsistent_datasets" => Dict{String,Any}(
            data.label => [s.deviation, s.uncertainty] for data in result.datasets for
            s in (_pair_sum(data.record),) if s !== nothing && s.consistent === false
        ),
        "flagged_mass_yields" => Dict{String,Any}(
            label => record.qualifiers for (label, record) in result.mass_yield_records if
            !isempty(_flagged(record.qualifiers))
        ),
        "retrievals" => _retrieval_runs(configuration),
    )
    level_density.model == "GC" && (
        record["shell_correction_file"] = _input_value(
            level_density.shell_correction_file, configuration
        )
    )
    configuration.yield_directory === nothing || (
        record["yield_directory"] = _relative_input(
            configuration.yield_directory, configuration
        )
    )
    configuration.yield_file === nothing ||
        (record["mass_yield_file"] = _relative_input(configuration.yield_file, configuration))
    return record
end

# The identities at the symmetric split; a test that has nothing to apply to is absent.
function _symmetry_record(symmetry::SymmetryDiagnostics)
    record = Dict{String,Any}(
        "pinned_curves" => symmetry.pinned_curves, "warnings" => symmetry.warnings
    )
    ismissing(symmetry.charge_set_invariant) ||
        (record["charge_set_invariant"] = symmetry.charge_set_invariant)
    ismissing(symmetry.R_T_at_symmetric_split) ||
        (record["R_T_at_symmetric_split"] = symmetry.R_T_at_symmetric_split)
    return record
end

# The correlogram, with NaN at a lag that holds no value, whether an autocorrelation could be
# fitted to it and, where it could, the coefficient the covariance was formed with, the scale of
# the trend's covariance against the expectation of χ², and the leave-one-out spread of ⟨R_T⟩
# over each yield distribution it was formed for.
function _trend_uncertainty_record(result::ExtractionResult)
    ρ = deviation_autocorrelation(result)
    record = Dict{String,Any}(
        "deviation_correlogram" =>
            Float64[coalesce(c, NaN) for c in result.deviation_correlogram],
        "autocorrelation_lags" => result.configuration.segments.autocorrelation_lags,
        "autocorrelation_estimated" => !ismissing(ρ),
    )
    ismissing(ρ) || (record["deviation_autocorrelation"] = ρ)
    index = findfirst(_is_trend, result.segmented_curves)
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

# One experiment pooled as one: its datasets, their accessions, the relation every record states
# where they agree, which decides how its curve was formed, and what each record states.
function _correlation_group_record(result::ExtractionResult, group::Vector{String})
    position = Dict(data.label => index for (index, data) in enumerate(result.datasets))
    indices = Int[position[label] for label in group]
    record = Dict{String,Any}(
        "datasets" => group,
        "accessions" => [_accession(result.datasets[i]) for i in indices],
        "relations_stated" =>
            [_correlation_relation(result.datasets[i].record) for i in indices],
    )
    relation = _relation(result.datasets, indices)
    isempty(relation) || (record["relation"] = relation)
    return record
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
# produced the files, from which of its configurations, and when. `load_configuration` has held
# every record to the version floor, so each carries one.
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
            runs[_relative_input(record, configuration)] = Dict{String,Any}(
                get(document, "run", Dict{String,Any}())
            )
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
    # A system with no recorded energy standard has no offset from one.
    offset.standard === nothing && return record
    record["standard_MeV"] = offset.standard
    record["standard_uncertainty_MeV"] = offset.standard_uncertainty
    record["offset_MeV"] = offset.offset
    return record
end

_parent(::Nothing) = nothing
_parent(path::AbstractString) = dirname(path)

# The ⟨TKE⟩(A) dataset that weighted the inversion: the file, its retrieval record entry, the
# parser revision that wrote it, a digest of the whole record, and the heavy masses that rest on
# a neighbour's value rather than their own; without one, why the inversion took no excitation
# weight.
function _kinetic_energy_record(result::ExtractionResult)
    configuration = result.configuration
    energies = result.mean_kinetic_energy
    weighted = manifest_domain(result).excitation_weighted
    energies === nothing && return Dict{String,Any}(
        "excitation_weighted" => weighted,
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
