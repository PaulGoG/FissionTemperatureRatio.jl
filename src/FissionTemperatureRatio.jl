"""
    FissionTemperatureRatio

Extraction of the temperature ratio `R_T = T_L/T_H` of complementary fully accelerated fission
fragments from experimental prompt neutron multiplicity data, described by joined straight
segments.

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
configuration = load_configuration("config/U233_nth.toml"; data_directory = "data")
result = run_pipeline(configuration)
write_results(result, "data/sims/U233_nth/" * run_identifier(configuration))
```

`scripts/run.jl` does this and writes the provenance record and the figures beside the tables.
"""
module FissionTemperatureRatio

using CSV: CSV
using DataFrames: DataFrame
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
using Measurements: Measurements
using SHA: sha1
using Statistics: mean, median, std
using TOML: TOML

include("retrieval.jl")
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

# What a run is driven by, reads and returns, exported.
export Configuration, load_configuration
export Multiplicity,
    read_multiplicity, read_multiplicity_directory, RatioCurve, multiplicity_ratio
export MeanKineticEnergy, read_mean_kinetic_energy
export RetrievalRecord, retrieval_record, retrieval_records
export SegmentedFit, fit_segments, segment_count, InsufficientDataError
export ExtractedCurve, TotalAverage, range_mean
export ExtractionResult, LeaveOneOut, run_pipeline, write_results, systematic_trend
export leave_one_out_spread, run_identifier, run_metadata
export publication_theme, plot_multiplicities, plot_ratio, save_figure, write_figures

# The rest of the interface, public and reached qualified: the builders and accessors a caller
# composing its own run needs, the stages of the method, and the constants they are defined by.
# Generic names such as `evaluate` and `covariance` are not exported, so that `using` does not
# bring them into a namespace where another package defines its own.
public FragmentationSettings, LevelDensitySettings, SegmentSettings, OutputSettings
public build_mass_table, build_level_density_model, build_charge_model
public A_H_range, has_symmetric_split, SHIPPED_MASS_TABLE, SHIPPED_SHELL_CORRECTIONS
public multiplicity, read_mass_yield_directory, mass_yield_coverage, yield_fraction
public mean_kinetic_energy_offset, PairSumScale, check_retrieval_versions, FLAGGED_QUALIFIERS
public DatasetDiagnostics, diagnose, consensus
public fit_weights, evaluate, covariance, pivots, unresolved_minimum
public CORRELOGRAM_LAGS, deviation_correlogram, autocorrelation_decay, pooled_correlation
public SymmetryDiagnostics, manifest_domain, curve_accessions, pooled_datasets
public superseded_datasets, deviation_autocorrelation, run_parameters

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exception, _, _
        stubs = (publication_theme, plot_multiplicities, plot_ratio, save_figure, write_figures)
        if exception.f in stubs &&
            Base.get_extension(@__MODULE__, :FissionTemperatureRatioCairoMakieExt) === nothing
            print(io, "\nThe figure functions are supplied by the CairoMakie extension: ")
            print(io, "`using CairoMakie` first.")
        end
        return nothing
    end
    return nothing
end

end
