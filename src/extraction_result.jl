# What a run returns: the result of an extraction and the records it holds.

"""
    SymmetryDiagnostics

Checks of the identities that hold by construction at the symmetric split.

# Fields

- `charge_set_invariant`: whether the charge numbers retained at the symmetric split are
  invariant under `Z -> Z₀ - Z`, which is what makes `r_ν = 1/2` invert to `R_T = 1` there;
  `missing` for a fissioning nucleus of odd mass number, which has no symmetric split, or a range
  that does not reach it.
- `R_T_at_symmetric_split`: the temperature ratio the inversion returns for `r_ν = 1/2` at
  `A₀/2`, one by identity for the charge-resolved inversion and for the ratio of means; `missing`
  where there is no symmetric split or the relation is undefined.
- `pinned_curves`: the labels of the segmented curves pinned to `r_ν = 1/2` at the symmetric
  split.
- `warnings`: every departure found, one sentence each; empty when the identities hold.
"""
struct SymmetryDiagnostics
    charge_set_invariant::Union{Bool,Missing}
    R_T_at_symmetric_split::Union{Float64,Missing}
    pinned_curves::Vector{String}
    warnings::Vector{String}
end

"""
    LeaveOneOut

The systematic trend refitted with one pooled dataset left out, or with the datasets of one
experiment left out together. The covariance of the trend follows the errors of a dataset along
the mass axis but, being built from deviations less their mean, holds no constant offset between
datasets; how far `⟨R_T⟩`, the number of segments and the breakpoints move when one of them is
removed is the measure of that.

# Fields

- `datasets`: the label of the pooled dataset left out, or the labels of the datasets of one
  experiment, which are pooled as one and left out as one.
- `outcome`: `"segmented curve"`, `"segmented curve, fitted without the required windows"`, or
  `"no fit: <reason>"`.
- `segments`, `breakpoints`: of the trend refitted without them; `0` and empty where there is no
  fit.
- `reduced_chi_squared`: `wrss/dof` of the refit; `NaN` where there is no fit.
- `total_average_R_T`: `⟨R_T⟩` of the refit over each yield distribution averaged over, keyed by
  the distribution's label; empty where there is no fit or no distribution.
"""
struct LeaveOneOut
    datasets::Vector{String}
    outcome::String
    segments::Int
    breakpoints::Vector{Int}
    reduced_chi_squared::Float64
    total_average_R_T::Dict{String,Float64}
end

"""
    JackknifeDesign

The refits of the systematic trend its jackknife uncertainty is formed over: the trend refitted
without each of a set of subsets of `d = k ÷ 2` of the `k` pooled datasets or experiments. With
two or three of them, `d = 1` and the refits are the [`LeaveOneOut`](@ref) refits.

# Fields

- `units`: `k`, the pooled datasets or experiments.
- `deleted`: `d = k ÷ 2`, the number left out of each refit.
- `subsets`: the refits asked for: every subset of `d` of the `k` where there are at most
  `segments.jackknife_subsets`, and that many drawn with `segments.jackknife_seed` otherwise.
- `complete`: whether those are every subset.
- `refits`: how many of them gave a curve; the jackknife is formed over these.
"""
struct JackknifeDesign
    units::Int
    deleted::Int
    subsets::Int
    complete::Bool
    refits::Int
end

"""
    ExtractionResult

Everything a run produces, held together so that it can be inspected interactively as well as
written to disk by [`write_results`](@ref).

# Fields

- `configuration`: the configuration the run was driven by.
- `datasets`: the experimental multiplicity data that was read.
- `r_ν`, `R_T`: the multiplicity and temperature ratios extracted point by point from each
  dataset.
- `masses`, `model`, `charge`, `domain`: the mass table, the level density model, the charge
  distribution and the fragmentation domain built from them, all FissionFragmentsDomain's.
- `averaging`: how the relation between `R_T` and `E*_H/TXE` was reduced over the charge
  distribution, with the excitation weights where a `⟨TKE⟩(A)` dataset was given.
- `mean_kinetic_energy`: that dataset carried onto the heavy-mass range, or `nothing`.
- `segmented_curves`: one [`ExtractedCurve`](@ref) per dataset that supports a fit and reaches the
  coverage floor, followed by the systematic-trend curve. A consuming code takes one of them.
- `range_mean_R_T`: for each segmented curve, the mean of its temperature ratio over its mass
  range with the uncertainty propagated through the curve's covariance, and for a trend that
  carries a jackknife the larger of that and the jackknife; see [`range_mean`](@ref) for why
  this is not the total average.
- `mass_yields`: the fragment mass yield distributions read, symmetrized where the configuration
  says so, empty when it names none.
- `mass_yield_records`: the [`RetrievalRecord`](@ref) of each of them, by label, where one lists
  it.
- `mass_yield_coverage`: for each of them, by label, the share of the heavy-fragment yield of the
  primary distribution at the masses it holds, [`mass_yield_coverage`](@ref). A distribution
  below `min_yield_coverage`, or excluded by the configuration, is read and reported but not
  averaged over.
- `total_average_R_T`: the quantity the literature quotes, `⟨R_T⟩` over each yield distribution
  that reaches the coverage floor, keyed by curve label and then by yield label. Empty when no
  yields were given.
- `consensus_r_ν`: the combined multiplicity ratio the systematic-trend curve was fitted to, so
  that the trend can be compared with its input.
- `dataset_diagnostics`: one record per multiplicity dataset read, in the order read.
- `dataset_outcomes`: for every dataset, whether it offers a segmented curve and, if not, why.
- `symmetry`: the identities at the symmetric split, checked.
- `segment_count_sensitivity`: for each yield distribution, `(segments, ⟨R_T⟩)` of the
  systematic trend at the selected number of segments and at one and two more, where the
  constraints admit them.
- `deviation_correlogram`: [`deviation_correlogram`](@ref) of the pooled datasets about the
  combined curve, lags one to [`CORRELOGRAM_LAGS`](@ref).
- `autocorrelation`: the coefficient fitted to its first lags, with which the covariance of the
  systematic trend was formed, [`deviation_autocorrelation`](@ref); `missing` where no lag could
  be estimated.
- `superseded`: the datasets a republication supersedes, by label, each with the label of its
  successor, through which alone it enters the pool; [`superseded_datasets`](@ref).
- `correlation_groups`: the labels of the pooled datasets that belong to one experiment, as
  their retrieval records name one another, one list per experiment. Such datasets are combined
  into one before they are pooled, so that the experiment counts once.
- `leave_one_out_requested`: whether the run was asked for the leave-one-out refits and the
  refits behind the jackknife of the systematic trend, [`run_pipeline`](@ref)'s `leave_one_out`.
- `leave_one_out`: one [`LeaveOneOut`](@ref) per pooled dataset or experiment, in the order read;
  empty with fewer than two of them or where the run was asked to skip the refits.
- `jackknife`: the [`JackknifeDesign`](@ref) of the jackknife the systematic trend carries;
  `nothing` where it carries none: with fewer than two pooled datasets or experiments, fewer than
  two refits that gave a curve, or a run asked to skip the refits.
- `curve_flags`: for each dataset curve that resolves no minimum over a range covering a
  required window, by label, the reason
  [`unresolved_minimum`](@ref) gives; such a curve stays in the manifest.
- `retrievals`: the retrieval records of the input directories as the run read them, by
  directory; the run record is written from these.
"""
struct ExtractionResult
    configuration::Configuration
    datasets::Vector{Multiplicity}
    r_ν::Vector{RatioCurve}
    R_T::Vector{RatioCurve}
    masses::MassExcessTable
    model::LevelDensityModel
    charge::ChargeModel
    domain::FragmentationDomain
    averaging::RatioAveraging
    mean_kinetic_energy::Union{MeanKineticEnergy,Nothing}
    segmented_curves::Vector{ExtractedCurve}
    range_mean_R_T::Dict{String,Tuple{Float64,Float64}}
    mass_yields::Vector{MassYield}
    mass_yield_records::Dict{String,RetrievalRecord}
    mass_yield_coverage::Dict{String,Float64}
    total_average_R_T::Dict{String,Dict{String,TotalAverage}}
    consensus_r_ν::RatioCurve
    dataset_diagnostics::Vector{DatasetDiagnostics}
    dataset_outcomes::Dict{String,String}
    symmetry::SymmetryDiagnostics
    segment_count_sensitivity::Dict{String,Vector{Tuple{Int,Float64}}}
    deviation_correlogram::Vector{Union{Missing,Float64}}
    autocorrelation::Union{Missing,Float64}
    superseded::Dict{String,String}
    correlation_groups::Vector{Vector{String}}
    leave_one_out_requested::Bool
    leave_one_out::Vector{LeaveOneOut}
    jackknife::Union{Nothing,JackknifeDesign}
    curve_flags::Dict{String,String}
    retrievals::RetrievalDocuments
end

"""
    manifest_domain(result) -> ManifestDomain

The `[domain]` record of a run: the level density model and its Gilbert-Cameron branch, the
ratio averaging and its excitation weighting, the charges per mass, the charge model, the mass
table, whether `ΔZ = 0` was imposed at the symmetric split and the version of
FissionFragmentsDomain, as that package spells them. A consuming code refuses a curve extracted
on a domain other than the one it partitions on.
"""
function manifest_domain(result::ExtractionResult)
    return ManifestDomain(
        result.model, result.averaging, result.domain, result.charge, result.masses
    )
end

"""
    systematic_trend(result) -> ExtractedCurve

The segmented curve fitted to the combination of the pooled datasets.
"""
function systematic_trend(result::ExtractionResult)
    index = findfirst(c -> c.label == SYSTEMATIC_TREND_LABEL, result.segmented_curves)
    index === nothing && throw(ArgumentError("this result carries no systematic-trend curve"))
    return result.segmented_curves[index]
end

# Why a curve offers no segmented description: the reason `fit_segments` gave, for a curve whose
# data cannot support a fit under the constraints, as for a dataset covering only a few mass
# pairs.
struct NoFit
    reason::String
end
