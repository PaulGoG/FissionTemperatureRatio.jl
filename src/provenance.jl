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
    "fragmentation.charge_distribution_file" => "dZfile",
    "fragmentation.fallback_charge_polarization" => "dZ",
    "fragmentation.fallback_charge_dispersion" => "sZ",
    "level_density.model" => "ldm",
    "level_density.ratio_averaging" => "avg",
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
required windows, the exclusion list, the charge distribution file, the yield directory — enters
as a content-hash token and is written in full into the run metadata, with the reason. The system
is not a token: it names the directory the identifier sits in. `significant_digits` changes how a
number is rendered, not the number, and is left out.
"""
function run_parameters(configuration::Configuration)
    fragmentation = configuration.fragmentation
    level_density = configuration.level_density
    segments = configuration.segments
    charge_file = fragmentation.charge_distribution_file
    yield_directory = configuration.yield_directory
    entries = Pair{String,Any}[
        "system.incident_energy" => configuration.system.incident_energy,
        "fragmentation.charges_per_mass" => fragmentation.charges_per_mass,
        "fragmentation.heavy_mass_min" => fragmentation.heavy_mass_min,
        "fragmentation.heavy_mass_max" => fragmentation.heavy_mass_max,
        "fragmentation.charge_distribution_file" => if charge_file === nothing
            "none"
        else
            _hash_token(_canonical(_relative_input(charge_file, configuration)))
        end,
        "fragmentation.fallback_charge_polarization" => fragmentation.fallback_charge_polarization,
        "fragmentation.fallback_charge_dispersion" => fragmentation.fallback_charge_dispersion,
        "level_density.model" => String(level_density.model),
        "level_density.ratio_averaging" => _averaging_name(level_density.ratio_averaging),
        "multiplicity.exclude" => _hash_token(_canonical(configuration.excluded_datasets)),
        "yield.subdirectory" => if yield_directory === nothing
            "none"
        else
            _hash_token(_canonical(_relative_input(yield_directory, configuration)))
        end,
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

# The values behind every hashed token, in full, for the run metadata.
function _hashed_values(configuration::Configuration)
    charge_file = configuration.fragmentation.charge_distribution_file
    yield_directory = configuration.yield_directory
    return Dict{String,Any}(
        RUN_IDENTIFIER_ABBREVIATIONS["fragmentation.charge_distribution_file"] =>
            charge_file === nothing ? "none" : _relative_input(charge_file, configuration),
        RUN_IDENTIFIER_ABBREVIATIONS["multiplicity.exclude"] =>
            Dict{String,Any}(configuration.excluded_datasets),
        RUN_IDENTIFIER_ABBREVIATIONS["yield.subdirectory"] => if yield_directory === nothing
            "none"
        else
            _relative_input(yield_directory, configuration)
        end,
        RUN_IDENTIFIER_ABBREVIATIONS["segments.required_windows"] =>
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

_averaging_name(::RatioOfMeans) = "ratio_of_means"
_averaging_name(::MeanOfRatios) = "mean_of_ratios"

# The system, as both the run metadata and the manifest state it. One definition, so the two
# cannot drift apart.
function _system_record(system::SystemSpecification)
    return Dict{String,Any}(
        "label" => system.label,
        "notation" => system_notation(system),
        "target_A" => system.target_A,
        "target_Z" => system.target_Z,
        "channel" => system.channel,
        "reaction" => system.reaction,
        "incident_energy_MeV" => system.incident_energy,
        "compound_A" => system.A₀,
        "compound_Z" => system.Z₀,
    )
end

"""
    run_metadata(result) -> Dict{String,Any}

What the library knows about a run, as nested tables ready for `TOML.print`: the system, every
configuration value in full, the identifier tokens and the values behind the hashed ones, the
input files, the package and dependency versions, and a summary of the result — each segmented
curve, the identities checked at the symmetric split, and what became of every dataset.

The caller adds what only it knows — the run identifier's place on disk, the configuration path,
the commit, the platform — and writes the whole; `scripts/run.jl` does.
"""
function run_metadata(result::ExtractionResult)
    configuration = result.configuration
    fragmentation = configuration.fragmentation
    level_density = configuration.level_density
    segments = configuration.segments
    relative(path) = path === nothing ? "none" : _relative_input(path, configuration)

    curves = Dict{String,Any}(
        curve.label => Dict{String,Any}(
            "segments" => FissionTemperatureRatio.segments(curve.fit),
            "breakpoints" => curve.fit.breakpoints,
            "pivots" => [[point[1], point[2]] for point in pivots(curve.fit)],
            "reduced_chi_squared" => curve.fit.wrss / curve.fit.dof,
            "bic" => curve.fit.bic,
            "bic_by_order" => [[order, value] for (order, value) in curve.fit.selection],
            "weights_imputed" => curve.fit.weights_imputed,
            "range_mean_R_T" => collect(result.range_mean_R_T[curve.label]),
        ) for curve in result.segmented_curves
    )
    symmetry = result.symmetry
    return Dict{String,Any}(
        "system" => _system_record(configuration.system),
        "configuration" => Dict{String,Any}(
            "charges_per_mass" => fragmentation.charges_per_mass,
            "heavy_mass_min" => fragmentation.heavy_mass_min,
            "heavy_mass_max" => fragmentation.heavy_mass_max,
            "fallback_charge_polarization" => fragmentation.fallback_charge_polarization,
            "fallback_charge_dispersion" => fragmentation.fallback_charge_dispersion,
            "model" => String(level_density.model),
            "ratio_averaging" => _averaging_name(level_density.ratio_averaging),
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
        "inputs" => Dict{String,Any}(
            "mass_excess_file" => relative(level_density.mass_excess_file),
            "shell_correction_file" => relative(level_density.shell_correction_file),
            "charge_distribution_file" => relative(fragmentation.charge_distribution_file),
            "multiplicity_directory" => relative(configuration.multiplicity_directory),
            "yield_directory" => relative(configuration.yield_directory),
            "datasets" => [basename(data.source) for data in result.datasets],
            "mass_yields" => [basename(data.source) for data in result.mass_yields],
        ),
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
                "R_a_at_symmetric_split" =>
                    if symmetry.R_a_at_symmetric_split === missing
                        "not applicable"
                    else
                        symmetry.R_a_at_symmetric_split
                    end,
                "pinned_curves" => symmetry.pinned_curves,
                "warnings" => symmetry.warnings,
            ),
        ),
    )
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
