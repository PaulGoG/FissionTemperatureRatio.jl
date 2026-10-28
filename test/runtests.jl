using TOML
using CSV: CSV
using DataFrames: DataFrame, nrow
using Test
using FissionFragmentsDomain
using Measurements: Measurements
using FissionTemperatureRatio
# The public names the tests reach that `using` does not bring in.
using FissionTemperatureRatio:
    A_H_range,
    has_symmetric_split,
    multiplicity,
    read_mass_yield_directory,
    mass_yield_coverage,
    yield_fraction,
    mean_kinetic_energy_offset,
    check_retrieval_versions,
    consensus,
    fit_weights,
    evaluate,
    covariance,
    pivots,
    unresolved_minimum,
    CORRELOGRAM_LAGS,
    deviation_correlogram,
    autocorrelation_decay,
    pooled_correlation,
    manifest_domain,
    curve_accessions,
    pooled_datasets,
    superseded_datasets,
    deviation_autocorrelation,
    run_parameters,
    RUN_IDENTIFIER_ABBREVIATIONS

include("fixtures.jl")

@testset "FissionTemperatureRatio" begin
    include("test_quality.jl")
    include("test_retrieval.jl")
    include("test_multiplicity_ratio.jl")
    include("test_kinetic_energy.jl")
    include("test_yields.jl")
    include("test_temperature_ratio.jl")
    include("test_segmented_fit.jl")
    include("test_consensus.jl")
    include("test_autocorrelation.jl")
    include("test_configuration.jl")
    include("test_pipeline.jl")
    include("test_experiments.jl")
    include("test_published_values.jl")
    include("test_manifest.jl")
    include("test_plotting.jl")
end
