# Entry point: julia --project scripts/run.jl config/<case>.toml

using DrWatson: projectdir
using FissionTemperatureRatio

function main(arguments::Vector{String})
    if length(arguments) != 1
        println("usage: julia --project scripts/run.jl <configuration.toml>")
        return 1
    end
    path = isabspath(arguments[1]) ? arguments[1] : joinpath(projectdir(), arguments[1])
    run_pipeline(load_configuration(path))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
