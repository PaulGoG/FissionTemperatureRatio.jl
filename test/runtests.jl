using TOML
using CSV: CSV
using DataFrames: DataFrame, nrow
using Test
using FissionFragmentsDomain
using Measurements: Measurements
using FissionTemperatureRatio

include("fixtures.jl")

@testset "FissionTemperatureRatio" begin
    include("test_quality.jl")
    include("test_multiplicity_ratio.jl")
    include("test_kinetic_energy.jl")
    include("test_yields.jl")
    include("test_temperature_ratio.jl")
    include("test_segmented_fit.jl")
    include("test_configuration.jl")
    include("test_pipeline.jl")
    include("test_published_values.jl")
    # Last: this one loads CairoMakie, and asserts the extension is absent until it does.
    include("test_manifest.jl")
    include("test_plotting.jl")
end
