# Run identification and metadata, so that every result traces back to a configuration, a commit
# and a machine.

"""
    run_identifier(configuration) -> String

A short identifier for a run, built from the parameters that change its result.

Two configurations that differ in any of these produce different identifiers, and two that agree
in all of them produce the same one, so results can be found again without consulting a log.
"""
function run_identifier(configuration::Configuration)
    parameters = Dict(
        "system" => configuration.system.label,
        "Z" => configuration.fragmentation.charges_per_mass,
        "AHmax" => configuration.fragmentation.A_H_max,
        "ldp" => String(configuration.level_density.prescription),
        "avg" => _averaging_name(configuration.level_density.ratio_averaging),
        "seg" => configuration.segments.max_segments,
        "minpts" => configuration.segments.min_points_per_segment,
        "pin" => configuration.segments.pin_symmetric_split,
    )
    return savename(parameters; connector = "_", sort = true)
end

_averaging_name(::RatioOfMeans) = "ratio_of_means"
_averaging_name(::MeanOfRatios) = "mean_of_ratios"

"""
    run_metadata(configuration) -> Dict{String,Any}

Metadata describing a run: the configuration it was driven by, the state of the source tree, and
the machine and Julia build that produced it.

Recorded alongside every set of results, so that a number can be attributed to a configuration, a
commit and a hardware platform without relying on memory.
"""
function run_metadata(configuration::Configuration)
    blas_threads = try
        LinearAlgebra.BLAS.get_num_threads()
    catch
        missing
    end
    cpu = Sys.cpu_info()

    return Dict{String,Any}(
        "run" => Dict{String,Any}(
            "identifier" => run_identifier(configuration),
            "timestamp" => Dates.format(Dates.now(), Dates.ISODateTimeFormat),
            "configuration_file" => configuration.source,
        ),
        "source" => Dict{String,Any}(
            "package_version" => string(PACKAGE_VERSION),
            "commit" => something(gitdescribe(projectdir()), "unavailable"),
        ),
        "julia" => Dict{String,Any}(
            "version" => string(VERSION),
            "threads" => Threads.nthreads(),
            "blas_threads" => blas_threads === missing ? "unavailable" : blas_threads,
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
