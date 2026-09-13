# Entry point: julia --project scripts/run.jl config/<system>.toml

using DrWatson: projectdir
using FissionTemperatureRatio

# Figures are supplied by the CairoMakie extension, a weak dependency, so that the package itself
# loads without a plotting stack. Load it where the environment provides it, and say plainly where
# it does not — the run is still complete, but it writes tables and metadata only.
if Base.identify_package("CairoMakie") === nothing
    @warn "CairoMakie is not reachable from this environment: the run will write its tables and \
           metadata but no figures. Install it into your default environment, which stays on the \
           load path alongside this project: `julia -e 'using Pkg; Pkg.add(\"CairoMakie\")'`."
else
    @eval using CairoMakie
end

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
