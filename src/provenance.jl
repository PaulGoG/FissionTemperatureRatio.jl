# Run identification and metadata, so that every result traces back to a configuration, a commit
# and a machine.

"""
    RUN_IDENTIFIER_ABBREVIATIONS

The token each configuration key is written as in a run identifier, keyed by the key's path in
the configuration file.

A run-identifier token is the configuration key it came from, so that a directory name can be
read back into a configuration without consulting the source. Spelling every key out in full
would put some two hundred characters into every run directory, so the keys are abbreviated —
and abbreviated **here**, in one exported table, rather than at the point of use. Every key
[`run_parameters`](@ref) uses has an entry, which the test suite asserts; a key without one is
an error rather than a silently invented token.
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
    "segments.max_segments" => "maxseg",
    "segments.min_points_per_segment" => "minpts",
    "segments.min_segment_span" => "minspan",
    "segments.pin_symmetric_split" => "pin",
    "segments.required_windows" => "win",
    "segments.windows_apply_to_datasets" => "windat",
    "segments.min_dataset_coverage" => "mincov",
)

# A value that does not reduce to a savename token — a list, a path — enters the identifier as a
# short content hash, and is written in full into the run metadata with this reason. The hash is
# of a canonical spelling, so the same content gives the same token on every machine.
const HASHED_TOKEN_REASON = "a list or a path does not reduce to a savename token; the token is \
                             the first eight hexadecimal digits of the SHA-1 of its canonical \
                             spelling, and the value is recorded here in full"

_hash_token(canonical::AbstractString) = bytes2hex(sha1(canonical))[1:8]

function _canonical(windows::Vector{UnitRange{Int}})
    return join(("$(first(w))-$(last(w))" for w in windows), ",")
end
_canonical(exclusions::Dict{String,String}) = join(sort!(collect(keys(exclusions))), "\n")
_canonical(path::AbstractString) = String(path)

# A path relative to the data directory, which is how the configuration spelt it.
function _relative_input(path::AbstractString, configuration::Configuration)
    return relpath(path, configuration.data_directory)
end

"""
    run_parameters(configuration) -> Dict{String,Any}

The configuration values that change the result, keyed by the token each contributes to the run
identifier.

Every key is abbreviated through [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref), and a key with no entry
there throws rather than being given a token on the spot. A value that is a list or a path — the
required windows, the exclusion list, a tabulated charge distribution, a mass or shell-correction
table other than the shipped one, the `⟨TKE⟩(A)` dataset, the yield directory — enters as a
content-hash token and is written in full into the run metadata, with the reason. The
Gilbert-Cameron branch and shell corrections are tokens of a Gilbert-Cameron run only, and read
`false` and `none` otherwise, as the run record states them. The system is not a token: it names
the directory the identifier sits in. `significant_digits` changes how a number is rendered, not
the number, and is left out.
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
        "fragmentation.charge_distribution_file" => _input_token(
            fragmentation.charge_distribution, configuration
        ),
        "fragmentation.zero_polarization_at_symmetry" => fragmentation.zero_polarization_at_symmetry,
        "level_density.model" => level_density.model,
        "level_density.mass_excess_file" => _input_token(
            level_density.mass_excess_file, configuration
        ),
        "level_density.shell_correction_file" => if gilbert_cameron
            _input_token(level_density.shell_correction_file, configuration)
        else
            "none"
        end,
        "level_density.deformed_branch" => gilbert_cameron && level_density.deformed_branch,
        "level_density.ratio_averaging" => level_density.ratio_averaging,
        "level_density.mean_kinetic_energy_file" => _input_token(
            level_density.mean_kinetic_energy_file, configuration
        ),
        "multiplicity.exclude" => _hash_token(_canonical(configuration.excluded_datasets)),
        "yield.subdirectory" => _input_token(configuration.yield_directory, configuration),
        "segments.max_segments" => segments.max_segments,
        "segments.min_points_per_segment" => segments.min_points_per_segment,
        "segments.min_segment_span" => segments.min_segment_span,
        "segments.pin_symmetric_split" => segments.pin_symmetric_split,
        "segments.required_windows" => _hash_token(_canonical(segments.required_windows)),
        "segments.windows_apply_to_datasets" => segments.windows_apply_to_datasets,
        "segments.min_dataset_coverage" => segments.min_dataset_coverage,
    ]
    parameters = Dict{String,Any}()
    for (key, value) in entries
        haskey(RUN_IDENTIFIER_ABBREVIATIONS, key) || throw(
            ArgumentError("no entry in RUN_IDENTIFIER_ABBREVIATIONS for the key $(repr(key))"),
        )
        parameters[RUN_IDENTIFIER_ABBREVIATIONS[key]] = value
    end
    return parameters
end

# A configured input: a named source ("wahl", "mean", "ame2020", "gilbert_cameron_1965") as itself, a path as the hash of
# its spelling relative to the data directory, an absent one as "none".
const NAMED_SOURCES = ("wahl", "mean", SHIPPED_MASS_TABLE, SHIPPED_SHELL_CORRECTIONS)

_input_token(::Nothing, ::Configuration) = "none"
function _input_token(value::AbstractString, configuration::Configuration)
    value in NAMED_SOURCES && return String(value)
    return _hash_token(_canonical(_relative_input(value, configuration)))
end

_input_value(::Nothing, ::Configuration) = "none"
function _input_value(value::AbstractString, configuration::Configuration)
    value in NAMED_SOURCES && return String(value)
    return _relative_input(value, configuration)
end

# The values behind every hashed token, in full, for the run metadata.
function _hashed_values(configuration::Configuration)
    level_density = configuration.level_density
    abbreviation(key) = RUN_IDENTIFIER_ABBREVIATIONS[key]
    return Dict{String,Any}(
        abbreviation("fragmentation.charge_distribution_file") =>
            _input_value(configuration.fragmentation.charge_distribution, configuration),
        abbreviation("level_density.mass_excess_file") =>
            _input_value(level_density.mass_excess_file, configuration),
        abbreviation("level_density.shell_correction_file") =>
            _input_value(level_density.shell_correction_file, configuration),
        abbreviation("level_density.mean_kinetic_energy_file") =>
            _input_value(level_density.mean_kinetic_energy_file, configuration),
        abbreviation("multiplicity.exclude") =>
            Dict{String,Any}(configuration.excluded_datasets),
        abbreviation("yield.subdirectory") =>
            _input_value(configuration.yield_directory, configuration),
        abbreviation("segments.required_windows") =>
            [[first(w), last(w)] for w in configuration.segments.required_windows],
    )
end

"""
    run_identifier(configuration) -> String

The token that names one run of one fissioning system: DrWatson's `savename` over
[`run_parameters`](@ref), tokens sorted and joined by `_`, every value rendered by one rule —
integers and integer-valued reals without a decimal part, everything else as printed.

Two configurations that differ in any result-changing key produce different identifiers, and two
that agree in all of them produce the same one, so a run can be found again without consulting a
log. The run directory is `data/sims/<system>/<identifier>/`.
"""
function run_identifier(configuration::Configuration)
    return savename(run_parameters(configuration); connector = "_", val_to_string = _token)
end

_token(x::Bool) = string(x)
_token(x::Integer) = string(x)
_token(x::AbstractString) = String(x)
# Trailing zeros carry no information and would make two spellings of one value.
_token(x::Real) = isinteger(x) ? string(Int(round(x))) : string(x)

"""
    run_metadata(result) -> Dict{String,Any}

What the library knows about a run, as nested tables ready for `TOML.print`: the system, every
configuration value in full, the identifier tokens and the values behind the hashed ones, the
fragmentation domain as the manifest records it together with the charge model that built it,
the input files and the retrievals that produced them — the parser revision and the record of the
`⟨TKE⟩(A)` dataset in particular, with the heavy masses it was interpolated or extrapolated to —
the package and dependency versions, and a summary of the result: each segmented curve, the
identities checked at the symmetric split, and what became of every dataset.

The caller adds what only it knows — the run identifier's place on disk, the configuration path,
the files written, the commit, the platform — and writes the whole; `scripts/run.jl` does.
"""
function run_metadata(result::ExtractionResult)
    configuration = result.configuration
    fragmentation = configuration.fragmentation
    level_density = configuration.level_density
    segments = configuration.segments
    relative(path) = path === nothing ? "none" : _relative_input(path, configuration)

    curves = Dict{String,Any}(
        curve.label => Dict{String,Any}(
            "kind" => curve.kind,
            "segments" => FissionTemperatureRatio.segments(curve.fit),
            "breakpoints" => curve.fit.breakpoints,
            "pivots" => [[point[1], point[2]] for point in pivots(curve.fit)],
            "reduced_chi_squared" => curve.fit.wrss / curve.fit.dof,
            "bic" => curve.fit.bic,
            "bic_by_order" => [[order, value] for (order, value) in curve.fit.selection],
            "weights_imputed" => curve.fit.weights_imputed,
            "pairs" => curve.pairs,
            "coverage" => curve.coverage,
            "range_mean_R_T" => collect(result.range_mean_R_T[curve.label]),
        ) for curve in result.segmented_curves
    )
    domain = manifest_domain(result)
    symmetry = result.symmetry
    flagged = Dict{String,Any}(
        label => tags for (label, tags) in result.qualifiers if !isempty(_flagged(tags))
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
            "max_segments" => segments.max_segments,
            "min_points_per_segment" => segments.min_points_per_segment,
            "min_segment_span" => segments.min_segment_span,
            "pin_symmetric_split" => segments.pin_symmetric_split,
            "required_windows" => [[first(w), last(w)] for w in segments.required_windows],
            "windows_apply_to_datasets" => segments.windows_apply_to_datasets,
            "min_dataset_coverage" => segments.min_dataset_coverage,
            "significant_digits" => configuration.output.significant_digits,
        ),
        "identifier" => Dict{String,Any}(
            "tokens" => run_parameters(configuration),
            "hashed" => _hashed_values(configuration),
            "hashed_reason" => HASHED_TOKEN_REASON,
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
            "zero_polarization_at_symmetry" => fragmentation.zero_polarization_at_symmetry,
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
            "datasets" => [basename(data.source) for data in result.datasets],
            "mass_yields" => [basename(data.source) for data in result.mass_yields],
            "flagged_datasets" => flagged,
            "flagged_mass_yields" => Dict{String,Any}(
                label => tags for
                (label, tags) in _yield_qualifiers(configuration, result.mass_yields) if
                !isempty(_flagged(tags))
            ),
            "retrievals" => _retrieval_runs(configuration),
        ),
        "mean_kinetic_energy" => _kinetic_energy_record(result),
        "source" => Dict{String,Any}(
            "package_version" => string(PACKAGE_VERSION),
            "dependencies" => _dependency_versions(),
        ),
        "result" => Dict{String,Any}(
            "segmented_curves" => curves,
            "dataset_outcomes" => Dict{String,Any}(result.dataset_outcomes),
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

# The `[run]` table of every retrieval record in the input directories: which parser revision
# produced the files, from which of its configurations, and when.
function _retrieval_runs(configuration::Configuration)
    runs = Dict{String,Any}()
    for directory in (
        configuration.multiplicity_directory,
        configuration.yield_directory,
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
    offset.standard === nothing && return Dict{String,Any}(
        "mean_MeV" => offset.mean, "standard" => "none recorded for this system"
    )
    return Dict{String,Any}(
        "mean_MeV" => offset.mean,
        "standard_MeV" => offset.standard,
        "standard_uncertainty_MeV" => offset.standard_uncertainty,
        "offset_MeV" => offset.offset,
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
    record = retrieval_record(energies.source)
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

# The manifests are not version-controlled, because this package supports a range of Julia
# versions and a manifest is resolved against one of them. The versions a run actually used are
# recorded here instead, so a result remains attributable to the code that produced it. The active
# project is read directly rather than through Pkg, which a library has no other reason to depend
# on.
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
