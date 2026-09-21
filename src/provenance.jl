# Run identification and metadata, so that every result traces back to a configuration, a commit
# and a machine.

"""
    RUN_IDENTIFIER_ABBREVIATIONS

The abbreviation each configuration key is written as in a run identifier.

A run-identifier token is the configuration key it came from, so that a file name can be read back
into a configuration without consulting the source. Spelling every key out in full would put some
hundred and eighty characters into every result file name, so the keys are abbreviated — and
abbreviated **here**, in one exported table, rather than at the point of use. Every key
[`run_identifier`](@ref) uses has an entry, which the test suite asserts; a key without one is an
error rather than a silently invented token.
"""
const RUN_IDENTIFIER_ABBREVIATIONS = Dict(
    "system" => "system",
    "charges_per_mass" => "nZ",
    "heavy_mass_min" => "AHmin",
    "heavy_mass_max" => "AHmax",
    "model" => "ldm",
    "ratio_averaging" => "avg",
    "max_segments" => "maxseg",
    "min_points_per_segment" => "minpts",
    "pin_symmetric_split" => "pin",
    "parsimony" => "parsimony",
    "incident_energy" => "E",
)

# The configuration keys whose values change the result, keyed as the configuration spells them.
# `required_windows` and the exclusion list are not here: neither reduces to a savename token, and
# both are recorded in full in the run metadata beside the result.
function _run_parameters(configuration::Configuration)
    parameters = Dict{String,Any}(
        "system" => configuration.system.label,
        "charges_per_mass" => configuration.fragmentation.charges_per_mass,
        "heavy_mass_min" => configuration.fragmentation.heavy_mass_min,
        "heavy_mass_max" => configuration.fragmentation.heavy_mass_max,
        "model" => String(configuration.level_density.model),
        "ratio_averaging" => _averaging_name(configuration.level_density.ratio_averaging),
        "max_segments" => configuration.segments.max_segments,
        "min_points_per_segment" => configuration.segments.min_points_per_segment,
        "pin_symmetric_split" => configuration.segments.pin_symmetric_split,
        "parsimony" => configuration.segments.parsimony,
    )
    # The label does not distinguish two incident energies of the same target and channel, which
    # are two systems; the identifier must.
    configuration.system.incident_energy > 0 &&
        (parameters["incident_energy"] = configuration.system.incident_energy)
    return parameters
end

"""
    run_identifier(configuration) -> String

A short identifier for a run, built from the configuration keys that change its result and written
with the abbreviations of [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref).

Two configurations that differ in any of those keys produce different identifiers, and two that
agree in all of them produce the same one, so results can be found again without consulting a log.
The required windows and the exclusion list are not among them; they are recorded in the run
metadata instead.

Throws an `ArgumentError` naming any key that has no abbreviation, rather than inventing one.
"""
function run_identifier(configuration::Configuration)
    parameters = Dict{String,Any}()
    for (key, value) in _run_parameters(configuration)
        haskey(RUN_IDENTIFIER_ABBREVIATIONS, key) || throw(
            ArgumentError("no entry in RUN_IDENTIFIER_ABBREVIATIONS for the key $(repr(key))"),
        )
        parameters[RUN_IDENTIFIER_ABBREVIATIONS[key]] = value
    end
    return savename(parameters; connector = "_", sort = true)
end

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
    run_metadata(configuration) -> Dict{String,Any}

Metadata describing a run: the configuration it was driven by, the state of the source tree, and
the machine and Julia build that produced it.

Recorded alongside every set of results, so that a number can be attributed to a configuration, a
commit and a hardware platform without relying on memory.
"""
function run_metadata(configuration::Configuration)
    blas_threads = LinearAlgebra.BLAS.get_num_threads()
    cpu = Sys.cpu_info()

    return Dict{String,Any}(
        "run" => Dict{String,Any}(
            "identifier" => run_identifier(configuration),
            "system" => _system_record(configuration.system),
            "segments" => Dict{String,Any}(
                "required_windows" => [
                    [first(w), last(w)] for w in configuration.segments.required_windows
                ],
                "windows_apply_to_datasets" =>
                    configuration.segments.windows_apply_to_datasets,
            ),
            "excluded_datasets" => Dict{String,Any}(configuration.excluded_datasets),
            "timestamp" => Dates.format(Dates.now(), Dates.ISODateTimeFormat),
            "configuration_file" => configuration.source,
        ),
        "source" => Dict{String,Any}(
            "package_version" => string(PACKAGE_VERSION),
            "dependencies" => _dependency_versions(),
            "commit" => something(gitdescribe(PACKAGE_ROOT), "unavailable"),
        ),
        "julia" => Dict{String,Any}(
            "version" => string(VERSION),
            "threads" => Threads.nthreads(),
            "blas_threads" => blas_threads,
            "versioninfo" => sprint(versioninfo),
        ),
        "platform" => Dict{String,Any}(
            "hostname" => gethostname(),
            "kernel" => string(Sys.KERNEL),
            "machine" => Sys.MACHINE,
            "cpu_model" => isempty(cpu) ? "unavailable" : cpu[1].model,
            "cpu_threads" => Sys.CPU_THREADS,
            "total_memory_GiB" => round(Sys.total_memory() / 2^30; digits = 2),
        ),
        "inputs" => Dict{String,Any}(
            "mass_excess_file" => configuration.level_density.mass_excess_file,
            "charge_distribution_file" => something(
                configuration.fragmentation.charge_distribution_file, "fallback only"
            ),
            "shell_correction_file" =>
                something(configuration.level_density.shell_correction_file, "not used"),
            "multiplicity_directory" => configuration.multiplicity_directory,
            "yield_directory" => something(configuration.yield_directory, "not used"),
        ),
    )
end

"""
    write_metadata(path, metadata)

Write run metadata to `path` as TOML, creating the directory if needed. Existing files are
preserved: a suffix is appended rather than overwriting a previous run's record.
"""
function write_metadata(path::AbstractString, metadata::AbstractDict)
    mkpath(dirname(path))
    target = _unused_path(path)
    open(target, "w") do io
        return TOML.print(io, metadata; sorted = true)
    end
    return target
end

# The resolved manifest of the active environment, copied beside the results. Manifests are not
# version-controlled, so this copy is what pins a run to exact dependency versions, the plotting
# stack included. Returns nothing when the active environment has no manifest.
function _write_environment_snapshot(directory::AbstractString, identifier::AbstractString)
    project = Base.active_project()
    project === nothing && return nothing
    candidates = (
        "Manifest-v$(VERSION.major).$(VERSION.minor).toml",
        "JuliaManifest.toml",
        "Manifest.toml",
    )
    index = findfirst(name -> isfile(joinpath(dirname(project), name)), candidates)
    index === nothing && return nothing
    target = _unused_path(joinpath(directory, "environment_$(identifier).toml"))
    mkpath(dirname(target))
    cp(joinpath(dirname(project), candidates[index]), target)
    return target
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
# rerun never destroys the record of the previous one.
function _unused_path(path::AbstractString)
    isfile(path) || return String(path)
    stem, extension = splitext(path)
    index = 1
    while isfile("$(stem)#$(index)$(extension)")
        index += 1
    end
    return "$(stem)#$(index)$(extension)"
end
