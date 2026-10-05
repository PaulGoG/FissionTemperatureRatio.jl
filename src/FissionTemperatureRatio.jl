"""
    FissionTemperatureRatio

Extraction of the temperature ratio `R_T = T_L/T_H` of complementary fully accelerated fission
fragments from experimental prompt neutron multiplicity data, and its parameterization by joined
straight segments.

Prompt emission model codes that share the total excitation energy at full acceleration take `R_T`
as input, either as one value for all fragmentations or as a function of the heavy-fragment mass
number. Obtaining it by fitting `ν(A)` with such a code makes the result dependent on the emission
model of that code. The method implemented here uses the same experimental `ν(A)` but no emission
calculation, so the ratio it produces is independent of any particular prompt emission treatment.

Two premises carry the extraction. The multiplicity ratio of complementary fragments follows their
excitation energy ratio, and the fragments are excited highly enough for their level densities to
be of Fermi-gas form, `E* = a T²`. Together these give

```
E*_L / E*_H = a_L T_L² / (a_H T_H²) ≈ ν_L / ν_H
```

so that for one fragmentation, with `ρ = a_L/a_H`, the heavy fragment carries
`E*_H/TXE = 1/(1 + ρ R_T²)`. The measured `r_ν = ν_H/(ν_L + ν_H)` is resolved by mass alone, and
`R_T(A_H)` is taken as the root of that relation reduced over the isobaric charge distribution at
`A_H`, every fragmentation with its own `ρ` and weighted by `p(Z, A_H)` and by its mean total
excitation. This is the exact inverse of a prompt emission code that partitions every fragment
pair with its own level density parameters. For a single effective ratio `R_a` it reduces to

```
R_T = [(1 - r_ν) / (R_a r_ν)]^(1/2),
```

eq. (4) of the paper, which remains selectable. The fragmentation domain, the level density
parameters and the relation itself are those of FissionFragmentsDomain.jl, shared with the codes
that consume `R_T(A_H)`.

The ratio extracted point by point from experimental data is scattered, and for some datasets
sparse, so it is the multiplicity ratio that is described by a continuous piecewise-linear
function whose segment count and breakpoints are selected from the data, and the temperature
ratio follows by the exact inversion above.

Method and conventions follow Eur. Phys. J. A **60**, 190 (2024). The vocabulary the package uses
for identifiers, configuration keys, files and column headers is set out in `docs/src/naming.md`.

# Entry point

```julia
configuration = load_configuration("config/U233_nth.toml")
result = run_pipeline(configuration)
write_results(result, "data/sims/U233_nth/" * run_identifier(configuration))
```

`scripts/run.jl` does this and writes the provenance record and the figures beside the tables.
"""
module FissionTemperatureRatio

using CSV: CSV
using DataFrames: DataFrame, eachrow
using Dates: Dates
using DrWatson: datadir, savename
using FissionFragmentsDomain:
    FissionFragmentsDomain,
    AME2020_MASS_EXCESS_FILE,
    ChargeDistribution,
    ChargeModel,
    ChargeResolved,
    FissioningSystem,
    FragmentationDomain,
    GILBERT_CAMERON_SHELL_CORRECTION_FILE,
    GilbertCameron,
    BackShiftedFermiGas,
    LevelDensityModel,
    MANIFEST_ABSCISSA,
    MANIFEST_ORDINATE,
    MEAN_KINETIC_ENERGY_SPEC,
    ManifestCurve,
    ManifestDomain,
    ManifestSystem,
    MassExcessTable,
    MassYield,
    Nuclide,
    RATIO_AVERAGINGS,
    RatioAveraging,
    SYSTEMATIC_TREND_LABEL,
    TemperatureRatioManifest,
    charge_model,
    charge_model_label,
    column,
    fragmentation_domain,
    integer_column,
    mass_yield,
    mean_charge_distribution,
    mean_total_excitation,
    neutron_induced_fission,
    read_charge_distribution,
    read_delimited_table,
    read_mass_excess_table,
    read_mass_yield,
    read_shell_correction_table,
    recommended_mean_total_kinetic_energy,
    spontaneous_fission,
    symmetric_charge_set_is_invariant,
    symmetrized_yield,
    system_record,
    temperature_ratio,
    temperature_ratio_slope,
    total_average,
    ratio_averaging,
    ratio_averaging_label,
    write_temperature_ratio_manifest
using LinearAlgebra: Symmetric, cond, dot, eigmin, tr
using Measurements: uncertainty, value
using SHA: sha1
using Statistics: mean, median, std
using TOML: TOML

# `@__DIR__` is resolved when this file is parsed, so the path is always available. `pkgdir` is
# not: it looks the module up in the loaded-package table, which is not yet populated while the
# module is still being defined, and returns `nothing` there. The project file is declared a
# dependency of the compiled module: the version is read when the module is compiled, and a
# release that changes nothing but the version would otherwise leave the earlier one in the cache
# and in the metadata of every run.
include_dependency(joinpath(dirname(@__DIR__), "Project.toml"))
const PACKAGE_VERSION = VersionNumber(
    TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))["version"]
)

include("multiplicity_ratio.jl")
include("kinetic_energy.jl")
include("yields.jl")
include("consensus.jl")
include("autocorrelation.jl")
include("temperature_ratio.jl")
include("segmented_fit.jl")
include("extracted_curve.jl")
include("configuration.jl")
include("plotting.jl")
include("pipeline.jl")
include("provenance.jl")

# Configuration
export Configuration, load_configuration
export build_mass_table, build_level_density_model, build_charge_model
export A_H_range, has_symmetric_split

# Input data
export Multiplicity, read_multiplicity, read_multiplicity_directory, multiplicity
export read_mass_yield_directory, mass_yield_coverage, yield_fraction
export MeanKineticEnergy, read_mean_kinetic_energy, mean_kinetic_energy_offset
export RetrievalRecord,
    retrieval_record,
    retrieval_qualifiers,
    check_retrieval_versions,
    FLAGGED_QUALIFIERS,
    pooling_weight,
    correlated_datasets,
    correlation_relation,
    pair_sum_scale

# Physics
export RatioCurve, multiplicity_ratio
export DatasetDiagnostics, diagnose, consensus

# Segmented description of the ratio
export SegmentedFit,
    fit_segments,
    fit_weights,
    evaluate,
    covariance,
    pivots,
    segments,
    unresolved_minimum,
    InsufficientDataError
export ExtractedCurve, TotalAverage, range_mean

# Pipeline
export ExtractionResult,
    SymmetryDiagnostics,
    run_pipeline,
    write_results,
    pool,
    systematic_trend,
    manifest_domain,
    curve_accessions,
    pooled_datasets,
    superseded_datasets,
    deviation_autocorrelation,
    CORRELOGRAM_LAGS,
    deviation_correlogram,
    autocorrelation_decay,
    pooled_correlation,
    LeaveOneOut,
    leave_one_out_spread
export run_identifier, run_parameters, run_metadata, RUN_IDENTIFIER_ABBREVIATIONS

# Figures
export publication_theme, plot_multiplicities, plot_ratio, save_figure, write_figures

end
