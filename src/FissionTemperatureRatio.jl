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

and hence, with `r_ν = ν_H/(ν_L + ν_H)` and `R_a = a_L/a_H`,

```
R_T = [(1 - r_ν) / (R_a r_ν)]^(1/2).
```

The ratio extracted point by point from experimental data is scattered, and for some datasets
sparse, so it is the multiplicity ratio that is described by a continuous piecewise-linear
function whose segment count and breakpoints are selected from the data, and the temperature
ratio follows by the exact transformation above.

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
using LinearAlgebra: Symmetric, cond, dot
using SHA: sha1
using Statistics: mean, median, std
using TOML: TOML

# `@__DIR__` is resolved when this file is parsed, so the path is always available. `pkgdir` is
# not: it looks the module up in the loaded-package table, which is not yet populated while the
# module is still being defined, and returns `nothing` there.
const PACKAGE_VERSION = VersionNumber(
    TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))["version"]
)

include("mass_data.jl")
include("level_density.jl")
include("fragmentation.jl")
include("multiplicity_ratio.jl")
include("yields.jl")
include("consensus.jl")
include("temperature_ratio.jl")
include("segmented_fit.jl")
include("segmented_curve.jl")
include("configuration.jl")
include("plotting.jl")
include("pipeline.jl")
include("provenance.jl")

# Configuration
export Configuration, load_configuration, build_level_density_model
export REACTIONS, CHANNEL_REACTION, system_label, system_notation, element_symbol
export A_H_range, has_symmetric_split

# Input data
export MassExcessTable, read_mass_excess_table, mass_excess
export ShellCorrectionTable, read_shell_correction_table
export ChargeDistribution, read_charge_distribution
export Multiplicity, read_multiplicity, read_multiplicity_directory, multiplicity
export MassYield, read_mass_yield, read_mass_yield_directory, mass_yield

# Physics
export LevelDensityModel, BackShiftedFermiGas, GilbertCameron
export level_density_parameter, shell_correction
export FragmentationDomain, fragmentation_domain, charges, charge_probability
export most_probable_charge, average_over_charge, symmetric_charge_set_is_invariant
export RatioAveraging, RatioOfMeans, MeanOfRatios
export RatioCurve, TREND_LABEL, multiplicity_ratio, level_density_ratio, temperature_ratio
export total_average
export DatasetDiagnostics, diagnose, consensus

# Segmented description of the ratio
export SegmentedFit,
    fit_segments, fit_weights, evaluate, covariance, pivots, segments, InsufficientDataError
export SegmentedCurve, TotalAverage, range_mean

# Pipeline
export ExtractionResult,
    SymmetryDiagnostics, run_pipeline, write_results, pool, systematic_trend
export run_identifier, run_parameters, run_metadata, RUN_IDENTIFIER_ABBREVIATIONS

# Figures
export publication_theme, plot_multiplicities, plot_ratio, save_figure, write_figures

end
