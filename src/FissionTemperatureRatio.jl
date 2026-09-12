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

The ratio extracted point by point from experimental data is scattered, and for some data sets
sparse, so it is the multiplicity ratio that is parameterized — by a continuous piecewise-linear
function whose segment count and breakpoints are selected from the data — and the temperature
ratio follows by the exact transformation above.

Method and conventions follow Eur. Phys. J. A **60**, 190 (2024).

# Entry point

```julia
configuration = load_configuration(joinpath(projectdir(), "config", "U233_nf.toml"))
result = run_pipeline(configuration)
```
"""
module FissionTemperatureRatio

using CSV: CSV
using CairoMakie
using DataFrames: DataFrame, eachrow
using Dates: Dates
using DrWatson: datadir, gitdescribe, projectdir, savename
using LaTeXStrings: @L_str
using LinearAlgebra: LinearAlgebra, Symmetric, cond, dot
using MathTeXEngine: texfont
using Statistics: median
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
include("temperature_ratio.jl")
include("segmented_fit.jl")
include("configuration.jl")
include("provenance.jl")
include("plotting.jl")
include("pipeline.jl")

# Configuration
export Configuration, load_configuration, build_prescription
export A_H_min, A_H_range, has_symmetric_split

# Input data
export MassExcessTable, read_mass_excess, mass_excess
export ShellCorrectionTable, read_shell_corrections
export ChargeDistributionData, read_charge_distribution
export MultiplicityData, read_multiplicity, read_multiplicity_directory, multiplicity

# Physics
export LevelDensityPrescription, BackShiftedFermiGas, GilbertCameron
export level_density_parameter, shell_correction
export FragmentationDomain, fragmentation_domain, charges, charge_probability
export most_probable_charge, average_over_charge, symmetric_charge_set_is_invariant
export RatioAveraging, RatioOfMeans, MeanOfRatios
export RatioCurve,
    TREND_LABEL, multiplicity_ratio, level_density_ratio, temperature_ratio, weighted_mean

# Parameterization
export SegmentedFit, fit_segments, fit_weights, evaluate, pivots, segments

# Pipeline
export PipelineResult, Parameterization, run_pipeline, write_results, pool, systematic_trend
export run_identifier, run_metadata

# Figures
export publication_theme, plot_multiplicities, plot_ratio, save_figure

end
