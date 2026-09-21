# Entry point: julia scripts/run.jl config/<system>.toml
#
# Runs in the environment of this directory, which carries CairoMakie, so a run always writes its
# figures and the plotting stack is resolved and recorded like every other dependency.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie: CairoMakie
using FissionTemperatureRatio

const ROOT = dirname(@__DIR__)

function main(arguments::Vector{String})
    if length(arguments) != 1
        println("usage: julia scripts/run.jl <configuration.toml>")
        return 1
    end
    # A path is taken as given when it exists from the working directory, and otherwise relative
    # to the repository root, so `config/<system>.toml` works from anywhere.
    path = isfile(arguments[1]) ? abspath(arguments[1]) : joinpath(ROOT, arguments[1])
    configuration = load_configuration(path; data_directory = joinpath(ROOT, "data"))
    run_pipeline(configuration; output_root = ROOT)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
