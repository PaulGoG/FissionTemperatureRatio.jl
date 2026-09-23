# Entry point: julia scripts/run.jl config/<system>.toml
#
# Runs in the environment of this directory, which carries CairoMakie, so a run always writes its
# figures and the plotting stack is resolved and recorded like every other dependency.
#
# Tables and the manifest go to data/sims/<system>/<run>/, figures to plots/<system>/<run>/, with
# the run named by its identifier. A run directory is never overwritten: one that already exists
# is moved aside as <run>#1, <run>#2, … before the new one is written. The provenance record
# metadata.toml sits beside the tables with copies of the configuration and of the resolved
# Manifest.toml of this environment. Set DRWATSON_STOREPATCH=true to record the uncommitted diff
# of a dirty checkout beside the commit.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie: CairoMakie
using FissionTemperatureRatio
using Dates: Dates
using DrWatson: tag!
using InteractiveUtils: versioninfo
using LinearAlgebra: BLAS
using Printf: @printf
using TOML: TOML

const ROOT = dirname(@__DIR__)

function main(arguments::Vector{String})
    if length(arguments) != 1
        println("usage: julia scripts/run.jl <configuration.toml>")
        return 1
    end
    # A path is taken as given when it exists from the working directory, and otherwise relative
    # to the repository root, so `config/<system>.toml` works from anywhere.
    path = isfile(arguments[1]) ? abspath(arguments[1]) : joinpath(ROOT, arguments[1])
    data = joinpath(ROOT, "data")
    configuration = load_configuration(path; data_directory = data)
    result = run_pipeline(configuration)

    system = configuration.system.label
    identifier = run_identifier(configuration)
    destination = joinpath(ROOT, "data", "sims", system, identifier)
    figures = joinpath(ROOT, "plots", system, identifier)
    move_aside(destination)
    move_aside(figures)
    write_results(result, destination)
    write_figures(result, figures)
    write_provenance(destination, result, path, data)

    for curve in result.segmented_curves
        averages = get(result.total_average_R_T, curve.label, nothing)
        averages === nothing && continue
        for (yield, entry) in sort!(collect(averages); by = first)
            @printf(
                "  ⟨R_T⟩ = %.4f ± %.4f  %s over Y(A) of %s\n",
                entry.value,
                entry.uncertainty,
                curve.label,
                yield
            )
        end
    end
    println("  written to ", destination)
    println("  figures in ", figures)
    return 0
end

# Never overwrite a run. An existing directory moves to <name>#1, <name>#2, …: the numbering
# DrWatson's safesave applies to a file, applied here to a whole run directory.
function move_aside(directory::AbstractString)
    isdir(directory) || return nothing
    index = 1
    while ispath("$(directory)#$(index)")
        index += 1
    end
    mv(directory, "$(directory)#$(index)")
    @info "existing run moved aside" from = directory to = "$(directory)#$(index)"
    return nothing
end

# The provenance record: what the library knows about the run, the identifier and the paths, the
# commit through DrWatson's tag!, and the platform, written as metadata.toml; beside it, copies of
# the configuration and of the resolved Manifest of this environment. A result then traces to
# configuration, commit, environment and hardware from its own directory.
function write_provenance(
    directory::AbstractString,
    result::ExtractionResult,
    configuration_path::AbstractString,
    data_directory::AbstractString,
)
    metadata = run_metadata(result)
    run = Dict{String,Any}(
        "identifier" => run_identifier(result.configuration),
        "configuration" => relpath(configuration_path, ROOT),
        "data_directory" => relpath(data_directory, ROOT),
        "script" => relpath(@__FILE__, ROOT),
        "written" => Dates.format(Dates.now(), Dates.ISODateTimeFormat),
    )
    tag!(run; gitpath = ROOT)
    metadata["run"] = run
    processors = Sys.cpu_info()
    metadata["platform"] = Dict{String,Any}(
        "julia" => string(VERSION),
        "hostname" => gethostname(),
        "kernel" => string(Sys.KERNEL),
        "machine" => Sys.MACHINE,
        "cpu" => isempty(processors) ? "unknown" : String(first(processors).model),
        "cpu_threads" => Sys.CPU_THREADS,
        "threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "total_memory_GiB" => round(Sys.total_memory() / 2^30; digits = 2),
        "versioninfo" => sprint(versioninfo),
    )
    open(joinpath(directory, "metadata.toml"), "w") do io
        return TOML.print(io, metadata; sorted = true)
    end
    cp(configuration_path, joinpath(directory, "configuration.toml"))
    manifest = joinpath(@__DIR__, "Manifest.toml")
    isfile(manifest) && cp(manifest, joinpath(directory, "Manifest.toml"))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
