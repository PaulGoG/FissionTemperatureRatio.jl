# The pipeline: from experimental multiplicity data to a segmented temperature ratio.

"""
    SymmetryDiagnostics

The identities that hold by construction at the symmetric split, checked rather than assumed.

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
    ExtractionResult

Everything a run produces, held together so that it can be inspected interactively as well as
written to disk by [`write_results`](@ref).

# Fields

- `configuration`: the configuration the run was driven by.
- `datasets`: the experimental multiplicity data that was read.
- `r_ν`, `R_T`: the multiplicity and temperature ratios extracted point by point from each
  dataset.
- `masses`, `model`, `charge`, `domain`: the mass table, the level density model, the charge
  distribution and the fragmentation domain built from them, all FissionFragmentsDomain's, which
  the consuming emission code partitions on too.
- `averaging`: how the relation between `R_T` and `E*_H/TXE` was reduced over the charge
  distribution, with the excitation weights where a `⟨TKE⟩(A)` dataset was given.
- `mean_kinetic_energy`: that dataset carried onto the heavy-mass range, or `nothing`.
- `qualifiers`: the reaction-code qualifiers the retrieval recorded for each multiplicity dataset,
  keyed by label.
- `segmented_curves`: one [`ExtractedCurve`](@ref) per dataset that supports a fit and reaches the
  coverage floor, followed by the systematic-trend curve. Alternatives offered to a prompt
  emission code, not an ensemble.
- `range_mean_R_T`: for each segmented curve, the mean of its temperature ratio over its mass
  range with the uncertainty propagated through the curve's covariance; see [`range_mean`](@ref)
  for why this is not the total average.
- `mass_yields`: the fragment mass yield distributions read, symmetrized where the configuration
  says so, empty when it names none.
- `mass_yield_coverage`: for each of them, by label, the share of the heavy-fragment yield of the
  primary distribution at the masses it holds, [`mass_yield_coverage`](@ref). A distribution
  below `min_yield_coverage`, or excluded by the configuration, is read and reported but not
  averaged over.
- `total_average_R_T`: the quantity the literature quotes, `⟨R_T⟩` over each yield distribution
  that reaches the coverage floor, keyed by curve label and then by yield label. Empty when no
  yields were given.
- `consensus_r_ν`: the combined multiplicity ratio the systematic-trend curve was fitted to, so
  that the trend can be checked against its own input rather than taken on trust.
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
- `correlation_groups`: the labels of the pooled datasets that belong to one experiment, as
  their retrieval records name one another, one list per experiment. Such datasets are combined
  into one before they are pooled, so that the experiment counts once.
- `leave_one_out`: one [`LeaveOneOut`](@ref) per pooled dataset or experiment, in the order read;
  empty with fewer than two of them or where the run was asked to skip them.
- `curve_flags`: for each dataset curve that resolves no minimum, by label, the reason
  [`unresolved_minimum`](@ref) gives; such a curve stays in the manifest.
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
    qualifiers::Dict{String,Vector{String}}
    segmented_curves::Vector{ExtractedCurve}
    range_mean_R_T::Dict{String,Tuple{Float64,Float64}}
    mass_yields::Vector{MassYield}
    mass_yield_coverage::Dict{String,Float64}
    total_average_R_T::Dict{String,Dict{String,TotalAverage}}
    consensus_r_ν::RatioCurve
    dataset_diagnostics::Vector{DatasetDiagnostics}
    dataset_outcomes::Dict{String,String}
    symmetry::SymmetryDiagnostics
    segment_count_sensitivity::Dict{String,Vector{Tuple{Int,Float64}}}
    deviation_correlogram::Vector{Union{Missing,Float64}}
    autocorrelation::Union{Missing,Float64}
    correlation_groups::Vector{Vector{String}}
    leave_one_out::Vector{LeaveOneOut}
    curve_flags::Dict{String,String}
end

"""
    manifest_domain(result) -> ManifestDomain

The `[domain]` record of a run: the level density model and its Gilbert-Cameron branch, the
ratio averaging and its excitation weighting, the charges per mass, the charge model, the mass
table and the version of FissionFragmentsDomain, as that package spells them. A consuming code
refuses a curve extracted on a domain other than the one it partitions on.
"""
function manifest_domain(result::ExtractionResult)
    return ManifestDomain(
        result.model, result.averaging, result.domain, result.charge, result.masses
    )
end

"""
    systematic_trend(result) -> ExtractedCurve

The segmented curve that follows the systematic behaviour of the ratio rather than any single
dataset.
"""
function systematic_trend(result::ExtractionResult)
    index = findfirst(c -> c.label == SYSTEMATIC_TREND_LABEL, result.segmented_curves)
    index === nothing && throw(ArgumentError("this result carries no systematic-trend curve"))
    return result.segmented_curves[index]
end

"""
    read_multiplicity_directory(directory) -> Vector{Multiplicity}

Read every `.dat` file in `directory`, sorted by name so that colours and markers are assigned
reproducibly. Other files are ignored, so a run record can sit beside the data it describes.

The label of each dataset is its file name stripped of the leading archive identifier and the
extension, with underscores replaced by spaces, which is how the measurements are named in the
literature. Two files of one author and year — two subentries of one measurement — would share
that label, and a label selects a curve downstream, so each of them carries its archive identifier
in parentheses instead.

Throws an `ArgumentError` where the directory does not exist, holds no data file, or holds two
files of one EXFOR accession: one measurement would be read, fitted and pooled twice.
"""
function read_multiplicity_directory(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("multiplicity directory not found: $(directory)"))
    files = _data_files(directory)
    isempty(files) &&
        throw(ArgumentError("multiplicity directory holds no data files: $(directory)"))
    _refuse_shared_accession(directory, files)
    return [
        read_multiplicity(joinpath(directory, file); label = label) for
        (file, label) in zip(files, _unique_labels(files))
    ]
end

# Data files carry the `.dat` extension; anything else in the directory — a run record, a note —
# is not data and is skipped rather than parsed and rejected.
function _data_files(directory::AbstractString)
    return sort!(filter!(f -> endswith(f, ".dat") && !startswith(f, "."), readdir(directory)))
end

# `<accession>_<Author>_<year>.dat`. The accession is eight digits, or nine where the dataset is a
# pointer within a subentry.
function _dataset_label(file::AbstractString)
    stem = replace(splitext(file)[1], r"^[0-9]+_" => "")
    label = replace(stem, '_' => ' ')
    # Initials are written without a space in the file names; restore it for the legend.
    return replace(label, r"(?<=\b[A-Z])\.(?=[A-Z][a-z])" => ". ")
end

function _accession(file::AbstractString)
    found = match(r"^([0-9]+)_", file)
    return found === nothing ? splitext(file)[1] : found.captures[1]
end

# The labels of the data files of one directory, the archive identifier in parentheses after each
# label two files share.
function _unique_labels(files::AbstractVector{<:AbstractString})
    labels = _dataset_label.(files)
    return [
        count(==(label), labels) == 1 ? label : "$(label) ($(_accession(file)))" for
        (file, label) in zip(files, labels)
    ]
end

# An EXFOR dataset identifier as the manifest admits it: entry, subentry and, where the subentry
# holds several datasets, a pointer character.
const ACCESSION_PATTERN = r"^[0-9A-Z]{8,9}$"

# The EXFOR accession of a dataset: the identifier its retrieval record gives, else the one its
# file name leads with, else none, for a tabulation that came from no archive.
function _dataset_accession(source::AbstractString)
    record = retrieval_record(source)
    candidate = if record !== nothing && haskey(record.entry, "identifier")
        string(record.entry["identifier"])
    else
        first(split(basename(source), '_'))
    end
    return occursin(ACCESSION_PATTERN, candidate) ? String(candidate) : ""
end

# Two files of one directory carrying one accession are one measurement held twice — a leftover
# of an earlier retrieval under another spelling of the author, a copy — and would be read, fitted
# and pooled twice, while an exclusion by that accession could name only one of them.
function _refuse_shared_accession(directory::AbstractString, files)
    held = Dict{String,String}()
    for file in files
        accession = _dataset_accession(joinpath(directory, file))
        isempty(accession) && continue
        haskey(held, accession) && throw(
            ArgumentError(
                "$(directory) holds two files of the accession $(accession), \
                 $(held[accession]) and $(file); one measurement would be read twice"
            ),
        )
        held[accession] = file
    end
    return nothing
end

# The exclusion key of every data file of a directory, with its label, from the file names and the
# retrieval record alone.
function _directory_keys(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("directory not found: $(directory)"))
    files = _data_files(directory)
    _refuse_shared_accession(directory, files)
    return Dict(
        _exclusion_key(joinpath(directory, file), label) => label for
        (file, label) in zip(files, _unique_labels(files))
    )
end

# The key a dataset is excluded under: its accession, or its label where it carries none.
function _exclusion_key(source::AbstractString, label::AbstractString)
    accession = _dataset_accession(source)
    return isempty(accession) ? String(label) : accession
end
_exclusion_key(data::Union{Multiplicity,MassYield}) = _exclusion_key(data.source, data.label)

# The exclusions that apply to `datasets`, by label, each with its reason. One that names no
# dataset read is an error: the configuration would claim an exclusion the run did not make.
function _excluded_labels(exclusions::Dict{String,String}, datasets, path::AbstractString)
    keys_read = Dict{String,String}()
    for data in datasets
        key = _exclusion_key(data)
        haskey(keys_read, key) && throw(
            ArgumentError(
                "$(repr(keys_read[key])) and $(repr(data.label)) are both held under \
                 $(repr(key)); an exclusion could not tell them apart"
            ),
        )
        keys_read[key] = data.label
    end
    for key in sort!(collect(keys(exclusions)))
        haskey(keys_read, key) || throw(
            ArgumentError(
                "$(path) excludes $(repr(key)), which names no dataset read; an exclusion \
                 is keyed by the EXFOR accession of its dataset"
            ),
        )
    end
    return Dict(keys_read[key] => reason for (key, reason) in exclusions)
end

"""
    pooled_datasets(result) -> Vector{String}

The labels of the datasets the systematic trend is combined from, in the order read: those that
provide at least one complete fragment pair within the fragmentation range and that the
configuration does not exclude. A dataset that forms no pair has nothing to pool, whether or not
an exclusion names it.
"""
function pooled_datasets(result::ExtractionResult)
    excluded = _excluded_labels(
        result.configuration.excluded_datasets, result.datasets, "multiplicity.exclude"
    )
    return [
        data.label for (data, curve) in zip(result.datasets, result.r_ν) if
        !isempty(curve) && !haskey(excluded, data.label)
    ]
end

"""
    deviation_correlogram(result) -> Vector{Union{Missing,Float64}}

The correlogram of the pooled datasets' deviations from the combined curve of a run, lags one to
[`CORRELOGRAM_LAGS`](@ref).
"""
deviation_correlogram(result::ExtractionResult) = result.deviation_correlogram

"""
    deviation_autocorrelation(result) -> Union{Float64,Missing}

The coefficient `ρ` with which the error of a pooled dataset is taken to be correlated along the
mass axis, `ρ^|A − A′|` between two of its mass numbers: [`autocorrelation_decay`](@ref) of
[`deviation_correlogram`](@ref) over the configured number of lags,
`segments.autocorrelation_lags`.

A dataset departs from the others by a slow drift along the mass axis rather than point by point;
with independent points the uncertainty of a smooth curve through the combined values would be
low by about `√((1 + ρ)/(1 − ρ))`. The covariance of the trend, and with it the uncertainty of
every tabulated `R_T` of the trend and of its total average, is formed with the correlation
matrix [`pooled_correlation`](@ref) builds from `ρ` and from the datasets combined at each mass
number. The curve of a single dataset is fitted to its points as independent ones.

`missing` where no lag could be estimated, the trend then being fitted with independent points.
"""
deviation_autocorrelation(result::ExtractionResult) = result.autocorrelation

"""
    leave_one_out_spread(result, mass_yield) -> Union{NamedTuple,Nothing}

The spread of the systematic trend's `⟨R_T⟩` over the yield distribution labelled `mass_yield`
across the [`LeaveOneOut`](@ref) refits of a run: `(; min, max, uncertainty)`, the least and the
greatest of the values `θᵢ`, `i = 1, …, k`, of the refits that hold one, and the delete-one
jackknife standard error over the pooled datasets,

```
σ = [(k − 1)/k Σᵢ (θᵢ − θ̄)²]^(1/2).
```

This is the uncertainty of the trend's `⟨R_T⟩` that comes from the datasets differing from one
another by more than their errors along the mass axis, which the covariance of the trend does not
contain. `nothing` where fewer than two refits hold a value.
"""
function leave_one_out_spread(result::ExtractionResult, mass_yield::AbstractString)
    θ = Float64[
        entry.total_average_R_T[mass_yield] for
        entry in result.leave_one_out if haskey(entry.total_average_R_T, mass_yield)
    ]
    k = length(θ)
    k < 2 && return nothing
    θ̄ = mean(θ)
    return (;
        min = minimum(θ), max = maximum(θ), uncertainty = sqrt((k - 1) / k * sum(abs2, θ .- θ̄))
    )
end

"""
    curve_accessions(result) -> Dict{String,String}

The EXFOR accession of every dataset of a run, keyed by its label: the identifier the retrieval
record states, or the one the file name leads with where no record lists the file, and empty for a
tabulation that carries neither. The systematic trend is fitted to several datasets and has
none; its entry is empty.

The label, `Author year`, is the display name and the key a consuming code selects a curve by; it
carries the accession only where two datasets would otherwise share it. The accession is what
identifies the measurement in the archive, so it is written into the manifest entry of every
dataset curve and into the name of every file of that dataset.
"""
function curve_accessions(result::ExtractionResult)
    accessions = Dict(data.label => _dataset_accession(data.source) for data in result.datasets)
    accessions[SYSTEMATIC_TREND_LABEL] = ""
    return accessions
end

# The token the tables and figures of each curve are named by, keyed by label. A dataset's is the
# stem of its input file, `<accession>_<Author>_<year>` for a retrieved one, so an output traces to
# its archive entry and to its input file by name alone; the trend's is its label.
function _file_tokens(result::ExtractionResult)
    tokens = Dict(
        data.label => _file_token(first(splitext(basename(data.source)))) for
        data in result.datasets
    )
    tokens[SYSTEMATIC_TREND_LABEL] = _file_token(SYSTEMATIC_TREND_LABEL)
    allunique(values(tokens)) || throw(
        ArgumentError(
            "two datasets would be written under one file name: $(sort!(collect(values(tokens))))",
        ),
    )
    return tokens
end

"""
    pool(curves; label) -> RatioCurve

Concatenate ratio curves from several datasets into one, sorted by mass number.

This is the naive combination, kept for comparison. It is **not** what builds the systematic-trend
curve: concatenating and weighting by the quoted uncertainties gives the result to whichever
author quoted the smallest ones and counts a dataset with many points more heavily than one with
few. [`consensus`](@ref) combines the datasets mass number by mass number instead, with a
between-dataset variance term.
"""
function pool(curves::Vector{RatioCurve}; label::AbstractString = "pooled")
    A_H = reduce(vcat, (curve.A_H for curve in curves); init = Int[])
    ratio = reduce(vcat, (curve.ratio for curve in curves); init = Float64[])
    σ = reduce(vcat, (curve.σ for curve in curves); init = Union{Missing,Float64}[])
    order = sortperm(A_H)
    return RatioCurve(A_H[order], ratio[order], σ[order], String(label))
end

"""
    run_pipeline(configuration; leave_one_out = true) -> ExtractionResult

Execute the extraction for one configuration.

The steps are: build the fragmentation domain from the mass table and the isobaric charge
distribution, and the level density model; read the experimental multiplicity data; form the
multiplicity ratio `r_ν = ν_H/(ν_L + ν_H)` of every dataset and the temperature ratio from it; and
describe the multiplicity ratio by joined straight segments, once per dataset and once more
following its systematic behaviour.

A run therefore returns several segmented curves rather than one. Where the datasets of a
fissioning nucleus disagree — as they do markedly for some — a curve fitted through all of them
describes none of them, and a prompt emission code has nothing to validate it against. They are
offered as alternatives: the code takes one as input, and which one describes reality is settled
by comparing the multiplicity distributions and yields it produces against experiment. The
systematic-trend curve is the one to use where a dataset is too sparse or too scattered to
determine the shape on its own.

The temperature ratio is obtained by inverting the fitted multiplicity ratio, not by fitting
segments to the temperature ratio a second time. The inversion is exact, whereas a second fit
would discard the uncertainty of the first and impose a piecewise-linear shape on a quantity that
is not piecewise-linear.

`leave_one_out` refits the systematic trend once with each pooled dataset left out,
[`LeaveOneOut`](@ref). The refits cost one trend fit per pooled dataset and can be skipped where
only the curves are wanted; the result then carries none.

Nothing is written. [`write_results`](@ref) writes the tables and the manifest of a result into a
directory of the caller's choosing; `scripts/run.jl` names that directory by the run identifier
and adds the provenance record and the figures.
"""
function run_pipeline(configuration::Configuration; leave_one_out::Bool = true)
    @info "reading input" configuration = configuration.source
    system = configuration.system
    settings = configuration.level_density
    masses = build_mass_table(settings)
    model = build_level_density_model(settings, masses)
    charge = build_charge_model(configuration, masses)
    range = A_H_range(configuration)
    domain = fragmentation_domain(
        system,
        charge,
        range;
        charges_per_mass = configuration.fragmentation.charges_per_mass,
        zero_polarization_at_symmetry = configuration.fragmentation.zero_polarization_at_symmetry,
    )
    @info "fragmentation domain" domain charge_model = charge_model_label(charge) level_density =
        settings.model

    mean_kinetic_energy = if settings.mean_kinetic_energy_file === nothing
        nothing
    else
        read_mean_kinetic_energy(settings.mean_kinetic_energy_file, system, range)
    end
    averaging = _averaging(settings.ratio_averaging, masses, domain, mean_kinetic_energy)

    datasets = read_multiplicity_directory(configuration.multiplicity_directory)
    # Before anything is fitted: an exclusion that names no dataset is a configuration the run
    # cannot honour. `load_configuration` has checked it; a configuration built by hand has not.
    excluded = _excluded_labels(
        configuration.excluded_datasets, datasets, "multiplicity.exclude"
    )
    if configuration.yield_directory === nothing
        # No directory, so the primary distribution alone is averaged over; an exclusion would be
        # written to the run record without having been applied.
        isempty(configuration.excluded_mass_yields) || throw(
            ArgumentError(
                "yield.exclude names \
                 $(join(sort!(collect(keys(configuration.excluded_mass_yields))), ", ")) \
                 without yield.subdirectory: the primary distribution alone is averaged over, \
                 and an exclusion has nothing to act on",
            ),
        )
    else
        held = _directory_keys(configuration.yield_directory)
        for key in sort!(collect(keys(configuration.excluded_mass_yields)))
            haskey(held, key) || throw(
                ArgumentError(
                    "yield.exclude excludes $(repr(key)), which names no distribution of \
                     $(configuration.yield_directory); an exclusion is keyed by the EXFOR \
                     accession of its dataset"
                ),
            )
        end
    end
    by_file = retrieval_qualifiers(configuration.multiplicity_directory)
    qualifiers = Dict(
        data.label => get(by_file, basename(data.source), String[]) for data in datasets
    )
    for data in datasets
        flagged = _flagged(qualifiers[data.label])
        isempty(flagged) ||
            @warn "the retrieval records a qualifier on this dataset; it is used, not \
                   corrected" dataset = data.label qualifiers = flagged
        scale = pair_sum_scale(retrieval_record(data.source))
        scale !== nothing &&
            scale.consistent === false &&
            @warn "the retrieval records the pair sum of this dataset off the scale of ν̄; a \
                   uniform scale cancels in r_ν, so it is used, not corrected" dataset =
                data.label deviation = scale.deviation uncertainty = scale.uncertainty
    end

    A₀ = system.compound.A
    r_ν = [multiplicity_ratio(data, A₀, range) for data in datasets]
    R_T = [temperature_ratio(averaging, model, domain, curve) for curve in r_ν]
    diagnostics = [diagnose(datasets[i], r_ν[i], A₀, range) for i in eachindex(datasets)]
    usable = findall(!isempty, r_ν)
    isempty(usable) &&
        throw(ArgumentError("no dataset provides both fragments of any pair within \
                       $(first(range)):$(last(range))"))

    segments_settings = configuration.segments
    segmented_curves = ExtractedCurve[]
    outcomes = Dict{String,String}()
    relation = (averaging, model, domain)

    # One curve per dataset. These are the alternatives a prompt emission code chooses between.
    # A dataset below the coverage floor is read, diagnosed and pooled, but offers no curve of
    # its own: the breakpoint search cannot place the minimum where the data do not reach.
    for (index, data) in enumerate(datasets)
        curve = r_ν[index]
        if isempty(curve)
            outcomes[data.label] = "no complete fragment pair within $(first(range)):$(last(range))"
            continue
        end
        coverage = diagnostics[index].coverage
        if coverage < segments_settings.min_pair_coverage
            outcomes[data.label] = "coverage $(round(coverage; digits = 3)) below the floor \
                                    $(segments_settings.min_pair_coverage); pooled only"
            @info "no segmented curve for this dataset" dataset = data.label coverage floor =
                segments_settings.min_pair_coverage
            continue
        end
        windows = if segments_settings.windows_apply_to_datasets
            segments_settings.required_windows
        else
            UnitRange{Int}[]
        end
        # Counted in measurements: an interpolated dataset's points count for raw over written.
        fraction = pooling_weight(retrieval_record(data.source))
        segmented = _segment(
            curve,
            relation,
            segments_settings,
            curve.label,
            windows,
            A₀;
            measured = fill(fraction, length(curve)),
        )
        if segmented isa ExtractedCurve
            push!(segmented_curves, segmented)
            outcomes[data.label] = "segmented curve"
        else
            outcomes[data.label] = "no fit: $(segmented)"
        end
    end
    unpinned = [c.label for c in segmented_curves if c.fit.pinned_value === nothing]
    if segments_settings.pin_symmetric_split && !isempty(unpinned)
        @info "fitted without the pin: no complete pair at the symmetric split" datasets =
            unpinned
    end
    # A dataset curve that resolves no minimum is offered all the same, and flagged.
    curve_flags = Dict{String,String}()
    for c in segmented_curves
        reason = unresolved_minimum(c.fit)
        reason === nothing && continue
        curve_flags[c.label] = reason
        @info "dataset curve flagged" dataset = c.label reason
    end

    # One curve following the systematic behaviour of the ratio: the minimum at the heavy magic
    # fragment placed rather than fitted, and the rise above the most probable fragmentation taken
    # through the whole body of data, so its slope falls between those of the individual datasets.
    # This is the curve to use where a dataset is too sparse or too scattered to determine the
    # shape on its own.
    # Datasets named in the configuration are kept out of the pooling but not out of the run: they
    # are still fitted, written and diagnosed, so an exclusion is visible rather than a silent
    # absence. An exclusion names its dataset by EXFOR accession.
    admitted = [i for i in usable if !haskey(excluded, datasets[i].label)]
    isempty(admitted) && throw(
        ArgumentError("every usable dataset is excluded from the pooling by configuration")
    )
    # Interpolated datasets count by their measured points, not their written rows.
    weights = [pooling_weight(retrieval_record(datasets[i].source)) for i in admitted]
    for (i, w) in zip(admitted, weights)
        w < 1 &&
            @info "interpolated dataset pooled at reduced weight" dataset = datasets[i].label weight =
                w
    end
    # Datasets of one experiment are combined into one before they are pooled.
    units = _pool_units(datasets, r_ν, admitted, weights)
    for unit in units.members
        length(unit) > 1 && @info "datasets of one experiment pooled as one" datasets = [
            datasets[i].label for i in unit
        ]
    end
    fitted = _fit_trend(units.curves, units.weights, relation, segments_settings, A₀)
    pooled = fitted.pooled
    combined = fitted.combined
    measured = fitted.measured
    windows = fitted.windows
    trend = fitted.trend
    trend isa ExtractedCurve && push!(segmented_curves, trend)

    isempty(segmented_curves) &&
        throw(ArgumentError("no dataset supports a segmented curve with \
                       min_points_per_segment = $(segments_settings.min_points_per_segment)"))

    for curve in segmented_curves
        @info "segmented curve" dataset = curve.label segments = segments(curve.fit) breakpoints =
            curve.fit.breakpoints reduced_chi_squared = curve.fit.wrss / curve.fit.dof
    end

    symmetry = _symmetry_diagnostics(configuration, relation, segmented_curves)
    for message in symmetry.warnings
        @warn message
    end

    # The primary distribution of the system: averaged over alone, or, beside a directory of
    # distributions, the reference their coverage is measured against.
    reference = if configuration.yield_file === nothing
        nothing
    else
        file = configuration.yield_file
        read_mass_yield(file; label = _dataset_label(basename(file)))
    end
    mass_yields = if configuration.yield_directory !== nothing
        read_mass_yield_directory(configuration.yield_directory)
    elseif reference !== nothing
        [reference]
    else
        MassYield[]
    end
    if configuration.symmetrize_yields && !isempty(mass_yields)
        # Before the qualifiers are looked up by file: the source and label are unchanged.
        mass_yields = [symmetrized_yield(y, A₀) for y in mass_yields]
        reference = symmetrized_yield(reference, A₀)
        @info "mass yields symmetrized: Y(A) and Y(A₀ - A) averaged where both are measured, \
               one wing standing for the other where it alone is"
    end
    for (distribution, tags) in _yield_qualifiers(configuration, mass_yields)
        flagged = _flagged(tags)
        isempty(flagged) ||
            @warn "the retrieval records a qualifier on this yield distribution; it is used, \
                   not corrected" distribution qualifiers = flagged
    end
    if mean_kinetic_energy !== nothing
        for distribution in mass_yields
            offset = mean_kinetic_energy_offset(mean_kinetic_energy, distribution, system)
            @info "⟨TKE⟩ of the excitation weights against the energy standard" dataset =
                mean_kinetic_energy.label yield = distribution.label mean_MeV = offset.mean standard_MeV =
                offset.standard offset_MeV = offset.offset
        end
    end
    # Coverage in yield, against the primary distribution: a set's own yields cannot measure it.
    # One below the floor, or excluded by the configuration, is reported but not averaged over.
    yield_coverage = Dict(
        y.label => mass_yield_coverage(y, reference, domain.heavy_masses) for y in mass_yields
    )
    excluded_yields = if configuration.yield_directory === nothing
        Dict{String,String}()
    else
        _excluded_labels(configuration.excluded_mass_yields, mass_yields, "yield.exclude")
    end
    averaged_yields = filter(mass_yields) do distribution
        coverage = yield_coverage[distribution.label]
        @info "yield coverage against the primary distribution" yield = distribution.label coverage reference =
            reference.label
        if haskey(excluded_yields, distribution.label)
            @info "no total average over this yield distribution: excluded by configuration" yield =
                distribution.label reason = excluded_yields[distribution.label]
            return false
        end
        coverage >= configuration.min_yield_coverage && return true
        @warn "no total average over this yield distribution" yield = distribution.label coverage floor =
            configuration.min_yield_coverage
        return false
    end
    total_average_R_T = _total_averages(segmented_curves, averaged_yields, domain.heavy_masses)
    for curve in segmented_curves
        for distribution in averaged_yields
            entry = get(
                get(total_average_R_T, curve.label, Dict()), distribution.label, nothing
            )
            entry === nothing && continue
            @info "total average" dataset = curve.label yield = distribution.label R_T =
                entry.value uncertainty = entry.uncertainty
        end
    end

    # The trend refitted with one and two segments more than selected: how far ⟨R_T⟩ rests on the
    # selected order.
    sensitivity = Dict{String,Vector{Tuple{Int,Float64}}}()
    if trend isa ExtractedCurve && !isempty(averaged_yields)
        selected = segments(trend.fit)
        alternatives = ExtractedCurve[trend]
        for extra in 1:2
            refit = _segment(
                combined,
                relation,
                segments_settings,
                SYSTEMATIC_TREND_LABEL,
                windows,
                A₀;
                measured = measured,
                order = selected + extra,
            )
            refit isa ExtractedCurve && push!(alternatives, refit)
        end
        for distribution in averaged_yields
            _overlaps(trend.R_T, distribution) || continue
            sensitivity[distribution.label] = [
                (segments(c.fit), first(total_average(c, distribution))) for c in alternatives
            ]
            @info "trend ⟨R_T⟩ against the number of segments" yield = distribution.label values = sensitivity[distribution.label]
        end
    end

    # The trend refitted with each pooled dataset left out in turn, the datasets of one experiment
    # together: how far ⟨R_T⟩, the number of segments and the breakpoints rest on one of them.
    # Only the values of a refit are used, so it takes the points as independent.
    refits = LeaveOneOut[]
    if leave_one_out && trend isa ExtractedCurve && length(units.curves) ≥ 2
        @info "leave-one-out refits of the systematic trend" refits = length(units.curves)
        for j in eachindex(units.curves)
            others = setdiff(eachindex(units.curves), j)
            refitted = _fit_trend(
                units.curves[others],
                units.weights[others],
                relation,
                segments_settings,
                A₀;
                estimate_autocorrelation = false,
                quiet = true,
            )
            left_out = [datasets[i].label for i in units.members[j]]
            refit = refitted.trend
            if !(refit isa ExtractedCurve)
                push!(
                    refits,
                    LeaveOneOut(
                        left_out, "no fit: $(refit)", 0, Int[], NaN, Dict{String,Float64}()
                    ),
                )
                continue
            end
            outcome = if refitted.windows == segments_settings.required_windows
                "segmented curve"
            else
                "segmented curve, fitted without the required windows"
            end
            averages = Dict{String,Float64}(
                distribution.label => first(total_average(refit, distribution)) for
                distribution in averaged_yields if _overlaps(refit.R_T, distribution)
            )
            push!(
                refits,
                LeaveOneOut(
                    left_out,
                    outcome,
                    segments(refit.fit),
                    copy(refit.fit.breakpoints),
                    refit.fit.wrss / refit.fit.dof,
                    averages,
                ),
            )
        end
    end

    return ExtractionResult(
        configuration,
        datasets,
        r_ν,
        R_T,
        masses,
        model,
        charge,
        domain,
        averaging,
        mean_kinetic_energy,
        qualifiers,
        segmented_curves,
        Dict(c.label => range_mean(c) for c in segmented_curves),
        mass_yields,
        yield_coverage,
        total_average_R_T,
        pooled,
        diagnostics,
        outcomes,
        symmetry,
        sensitivity,
        fitted.correlogram,
        fitted.autocorrelation,
        [[datasets[i].label for i in unit] for unit in units.members if length(unit) > 1],
        refits,
        curve_flags,
    )
end

# One field of the pair-sum scale the retrieval records for a dataset, or `missing`.
function _scale_field(data::Multiplicity, field::Symbol, digits::Integer)
    scale = pair_sum_scale(retrieval_record(data.source))
    scale === nothing && return missing
    value = getfield(scale, field)
    return value isa AbstractFloat ? round(value; sigdigits = digits) : value
end

# The reaction-code qualifiers of the yield distributions a run read, keyed by label; the record
# lists every dataset its query accepted, most of which need not be staged.
function _yield_qualifiers(configuration::Configuration, mass_yields::Vector{MassYield})
    isempty(mass_yields) && return Dict{String,Vector{String}}()
    by_file = Dict{String,Vector{String}}()
    for directory in unique(dirname(y.source) for y in mass_yields)
        merge!(by_file, retrieval_qualifiers(directory))
    end
    return Dict(
        distribution.label => get(by_file, basename(distribution.source), String[]) for
        distribution in mass_yields
    )
end

# The averaging a configuration names. The charge-resolved inversion is weighted by the mean total
# excitation of every fragmentation where a ⟨TKE⟩(A) dataset is given, and by p(Z, A_H) alone
# otherwise; the run says which.
function _averaging(
    label::AbstractString,
    masses::MassExcessTable,
    domain::FragmentationDomain,
    mean_kinetic_energy::Union{MeanKineticEnergy,Nothing},
)
    label == "charge_resolved" || return ratio_averaging(label)
    if mean_kinetic_energy === nothing
        @info "no ⟨TKE⟩(A) dataset: the charge-resolved inversion weights each fragmentation by \
               p(Z, A_H) alone, not by its excitation"
        return ChargeResolved()
    end
    record = retrieval_record(mean_kinetic_energy.source)
    flagged = record === nothing ? String[] : _flagged(record.qualifiers)
    isempty(flagged) ||
        @warn "the retrieval records a qualifier on the ⟨TKE⟩(A) dataset; it is used, not \
               corrected" dataset = mean_kinetic_energy.label qualifiers = flagged
    @info "the charge-resolved inversion is weighted by ⟨TXE⟩ = Q + E*_CN - ⟨TKE⟩(A_H)" dataset =
        mean_kinetic_energy.label interpolated = mean_kinetic_energy.interpolated extrapolated =
        mean_kinetic_energy.extrapolated
    return ChargeResolved(mean_total_excitation(masses, domain, mean_kinetic_energy.values))
end

# Whether a yield distribution carries positive weight at any mass number of the curve.
function _overlaps(curve::RatioCurve, distribution::MassYield)
    total = 0.0
    for mass in curve.A_H
        entry = mass_yield(distribution, mass)
        entry === nothing || (total += entry[1])
    end
    return total > 0
end

# The total average of every segmented curve over every yield distribution. A pair that shares no
# mass number is omitted rather than reported as zero: a distribution covering only the light wing
# says nothing about a ratio defined on the heavy one.
function _total_averages(
    curves::Vector{ExtractedCurve},
    mass_yields::Vector{MassYield},
    heavy_masses::AbstractUnitRange{<:Integer},
)
    averages = Dict{String,Dict{String,TotalAverage}}()
    isempty(mass_yields) && return averages
    for curve in curves
        per_distribution = Dict{String,TotalAverage}()
        for distribution in mass_yields
            if _overlaps(curve.R_T, distribution)
                per_distribution[distribution.label] = TotalAverage(
                    curve, distribution, heavy_masses
                )
            else
                @warn "no total average" dataset = curve.label yield = distribution.label reason = "no mass number in common carries a positive yield"
            end
        end
        isempty(per_distribution) || (averages[curve.label] = per_distribution)
    end
    return averages
end

# Fit one ratio curve and carry it through to the temperature ratio. Returns the reason, as a
# string, when the curve cannot support a fit at all, which happens for datasets covering only a
# few mass pairs; `quiet` keeps that from the log where the caller reports it otherwise.
#
# `fit_segments` pins at the first abscissa it is given. r_ν = 1/2 is an identity at A₀/2 only, so
# the pin is applied to a curve whose first complete pair is the symmetric split and to no other;
# a dataset starting above it is fitted unpinned.
function _segment(
    curve::RatioCurve,
    relation::Tuple{RatioAveraging,LevelDensityModel,FragmentationDomain},
    settings::SegmentSettings,
    label::AbstractString,
    windows::Vector{UnitRange{Int}},
    A₀::Integer;
    measured::AbstractVector{<:Real} = ones(length(curve)),
    order::Union{Integer,Nothing} = nothing,
    correlation::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    quiet::Bool = false,
)
    pinned = settings.pin_symmetric_split && !isempty(curve) && 2 * first(curve.A_H) == A₀
    fit = try
        fit_segments(
            curve.A_H,
            curve.ratio,
            curve.σ;
            max_segments = order === nothing ? settings.max_segments : order,
            min_segments = order === nothing ? 1 : order,
            measured = measured,
            min_points_per_segment = settings.min_points_per_segment,
            min_segment_span = settings.min_segment_span,
            pinned_value = pinned ? 0.5 : nothing,
            required_windows = windows,
            bounds = (0.0, 1.0),
            correlation = correlation,
        )
    catch exception
        exception isa InsufficientDataError || rethrow()
        order === nothing &&
            !quiet &&
            @warn "no segmented curve for this dataset" dataset = label reason =
                exception.msg
        return exception.msg
    end
    return ExtractedCurve(label, fit, curve, relation...)
end

# What a trend is combined from: each admitted dataset, or, for the datasets of one experiment —
# those whose retrieval records name one another under `correlated_with` — their combination,
# formed as a pool is and entering the pool as one curve at the uncertainty of one measurement
# with the measured fraction of each of its points. Runs or analyses on one apparatus share their
# systematic errors; pooled side by side they would count as independent measurements in the
# between-dataset variance and weigh as several. `members` holds, for each curve, the indices of
# its datasets.
function _pool_units(
    datasets::Vector{Multiplicity},
    r_ν::Vector{RatioCurve},
    admitted::Vector{Int},
    weights::Vector{Float64},
)
    accessions = [_dataset_accession(datasets[i].source) for i in admitted]
    partners = [correlated_datasets(retrieval_record(datasets[i].source)) for i in admitted]
    # The experiments: the connected components of the relation, by position in `admitted`.
    experiment = collect(eachindex(admitted))
    for a in eachindex(admitted), b in (a + 1):length(admitted)
        named =
            (!isempty(accessions[b]) && accessions[b] in partners[a]) ||
            (!isempty(accessions[a]) && accessions[a] in partners[b])
        named || continue
        merged, kept = experiment[b], experiment[a]
        replace!(experiment, merged => kept)
    end

    curves = RatioCurve[]
    unit_weights = Union{Float64,Vector{Float64}}[]
    members = Vector{Int}[]
    for first_position in unique(experiment)
        positions = findall(==(first_position), experiment)
        push!(members, admitted[positions])
        if length(positions) == 1
            push!(curves, r_ν[admitted[first_position]])
            push!(unit_weights, weights[first_position])
            continue
        end
        label = join((datasets[i].label for i in admitted[positions]), " + ")
        combined, measured, σ_measurement, _, _ = _consensus(
            r_ν[admitted[positions]], label, weights[positions]
        )
        push!(curves, RatioCurve(combined.A_H, combined.ratio, σ_measurement, label))
        push!(unit_weights, measured)
    end
    return (; curves, weights = unit_weights, members)
end

# The systematic trend of a pool of datasets: their combined curve, the correlogram of their
# deviations from it, the autocorrelation fitted to that, and the segmented curve, fitted with the
# required windows and, where no curve satisfies them, without. The trend is the curve, or the
# reason there is none; `windows` are those it was fitted with, and `autocorrelation` is
# `missing` where no lag of the correlogram could be estimated. Without `estimate_autocorrelation`
# the points are taken as independent, which suffices where only the values of the fit are used.
function _fit_trend(
    curves::Vector{RatioCurve},
    weights::AbstractVector,
    relation::Tuple{RatioAveraging,LevelDensityModel,FragmentationDomain},
    settings::SegmentSettings,
    A₀::Integer;
    estimate_autocorrelation::Bool = true,
    quiet::Bool = false,
)
    pooled, measured, σ_measurement, members, shares = _consensus(
        curves, SYSTEMATIC_TREND_LABEL, weights
    )
    # The trend is fitted to the combined values at the uncertainty of one measurement, with the
    # measured fraction beside it, as a dataset is: the standard error of `pooled` carries the
    # fraction already, and fitting to it would count the fraction twice.
    combined = RatioCurve(pooled.A_H, pooled.ratio, σ_measurement, pooled.label)
    # A curve that is itself a combination enters the correlogram with its mean factor.
    factors = Float64[w isa Real ? w : mean(w) for w in weights]
    correlogram = deviation_correlogram(curves, pooled; weights = factors)
    lags = min(settings.autocorrelation_lags, length(correlogram))
    estimated = estimate_autocorrelation && !all(ismissing, view(correlogram, 1:lags))
    ρ = estimated ? autocorrelation_decay(correlogram; lags = lags) : 0.0
    # The errors of the combined points: correlated along the mass axis within a dataset,
    # independent between datasets.
    correlation = ρ > 0 ? pooled_correlation(pooled.A_H, members, shares, ρ) : nothing
    windows = settings.required_windows
    trend = _segment(
        combined,
        relation,
        settings,
        SYSTEMATIC_TREND_LABEL,
        windows,
        A₀;
        measured = measured,
        correlation = correlation,
        quiet = quiet,
    )
    if !(trend isa ExtractedCurve) && !isempty(windows)
        quiet ||
            @warn "no segmented curve satisfies the required windows; the systematic-trend \
                   curve was fitted without them" windows reason = trend
        windows = UnitRange{Int}[]
        trend = _segment(
            combined,
            relation,
            settings,
            SYSTEMATIC_TREND_LABEL,
            windows,
            A₀;
            measured = measured,
            correlation = correlation,
            quiet = quiet,
        )
    end
    autocorrelation = estimated ? ρ : missing
    return (; pooled, combined, measured, correlogram, autocorrelation, trend, windows)
end

# The identities that hold by construction at the symmetric split, checked rather than assumed.
function _symmetry_diagnostics(
    configuration::Configuration,
    relation::Tuple{RatioAveraging,LevelDensityModel,FragmentationDomain},
    curves::Vector{ExtractedCurve},
)
    averaging, model, domain = relation
    warnings = String[]
    A₀ = configuration.system.compound.A
    invariant = something(symmetric_charge_set_is_invariant(domain), missing)
    if invariant === false
        push!(
            warnings,
            "the charge numbers retained at the symmetric split are not invariant under \
             Z -> Z₀ - Z, so r_ν = 1/2 does not return R_T = 1 there exactly; this occurs when \
             the charge distribution is polarized at A₀/2, or when Z₀ is odd",
        )
    end

    R_T_symmetric = missing
    if iseven(A₀) && (A₀ ÷ 2) in domain.heavy_masses
        value = temperature_ratio(averaging, model, domain, A₀ ÷ 2, 0.5)
        if value !== nothing
            R_T_symmetric = value
            deviation = abs(value - 1)
            deviation < 1e-8 || push!(
                warnings,
                "r_ν = 1/2 inverts to R_T = $(round(value; sigdigits = 8)) at the symmetric \
                 split under ratio_averaging = \"$(ratio_averaging_label(averaging))\", \
                 $(round(deviation; sigdigits = 3)) from unity where the two fragments are the \
                 same nuclide",
            )
        end
    end

    pinned = String[]
    for c in curves
        c.fit.pinned_value === nothing && continue
        push!(pinned, c.label)
        # A pin anywhere but the symmetric split asserts r_ν = 1/2 where no identity holds.
        2 * c.fit.x₀ == A₀ || push!(
            warnings,
            "the curve \"$(c.label)\" is pinned at A_H = $(c.fit.x₀), which is not the \
             symmetric split of A₀ = $(A₀)",
        )
    end
    return SymmetryDiagnostics(invariant, R_T_symmetric, pinned, warnings)
end

# The run record a consuming code reads: the system, the quantity tabulated, the domain the curves
# were extracted on, and for every segmented curve its label, its kind, the two files it is
# tabulated in and, for a dataset curve, its archive accession, written by FissionFragmentsDomain's
# writer and nothing besides. What else a run
# knows about each curve is in the segmented-curves table and in the run metadata.
function _write_manifest(
    result::ExtractionResult,
    directory::AbstractString,
    identifier::AbstractString,
    written::Dict{String,String},
)
    record = system_record(result.configuration.system)
    system = ManifestSystem(
        record["label"],
        record["notation"],
        record["target_A"],
        record["target_Z"],
        record["channel"],
        record["reaction"],
        record["incident_energy_MeV"],
        record["compound_A"],
        record["compound_Z"],
    )
    accessions = curve_accessions(result)
    curves = [
        ManifestCurve(
            curve.label,
            curve.kind,
            basename(written["R_T_segmented/$(curve.label)"]),
            basename(written["r_nu_pivots/$(curve.label)"]);
            accession = accessions[curve.label],
        ) for curve in result.segmented_curves
    ]
    path = joinpath(directory, "manifest_$(identifier).toml")
    manifest = TemperatureRatioManifest(
        system,
        MANIFEST_ORDINATE,
        copy(MANIFEST_ABSCISSA),
        # Stated so that a reader knows what the columns hold, not so that a reader looks them up
        # by name: the files are read by column position.
        ["A_H", "R_T", "R_T_uncertainty"],
        curves,
        path,
        manifest_domain(result),
    )
    return write_temperature_ratio_manifest(path, manifest)
end

function _summarize_range(masses::Vector{Int})
    isempty(masses) && return "none"
    return "$(first(masses)):$(last(masses)) ($(length(masses)) mass numbers)"
end

# A ratio table: the abscissa, the quantity, then its uncertainty. The quantity names its own
# column, so a reader holding the file knows what is in it without consulting the file name. An
# unquoted uncertainty is an empty field, never a zero, which would denote an exact value.
function _ratio_table(curve::RatioCurve, quantity::AbstractString, significant_digits::Integer)
    return DataFrame(
        :A_H => curve.A_H,
        Symbol(quantity) => round.(curve.ratio; sigdigits = significant_digits),
        Symbol("$(quantity)_uncertainty") => [
            ismissing(σ) ? missing : round(σ; sigdigits = significant_digits) for σ in curve.σ
        ],
    )
end

"""
    write_results(result, directory) -> Dict{String,String}

Write the tabulated ratios, the segment pivots, the segmented-curve summary, the total averages,
the leave-one-out refits of the systematic trend, the dataset diagnostics and the manifest of a
run into `directory`, and return the paths written, keyed by content.

`directory` is created. One that already holds files is refused rather than written into, so a
caller that wants a second run of one configuration beside the first moves the first aside;
`scripts/run.jl` does, numbering it `#1`, `#2`, …, and names the directory
`data/sims/<system>/<run identifier>/`. The provenance record and the figures are not written
here: [`run_metadata`](@ref) supplies what the library knows about a run, and the script adds the
rest, writes `metadata.toml`, and calls [`write_figures`](@ref).

The manifest, `manifest_<run identifier>.toml`, is the contract with a consuming code, which
stages the whole directory and selects the manifest by that prefix; it holds the system, the
domain and, per segmented curve, the label, the kind, the two files and the accession of a dataset
curve, and nothing else. What else
a run reports about each curve — segments, pin, span, pairs, coverage, reduced chi-squared, range
mean — is one row per manifest curve, keyed by its label, in `segmented_curves_<run
identifier>.csv`; the total averages are in `total_average_R_T_<run identifier>.csv`, and the
trend refitted with each pooled dataset left out, where the run made those refits, in
`leave_one_out_<run identifier>.csv`. `dataset_diagnostics.csv` flags, with the reason, every
dataset curve that resolves no minimum.

The tables of a dataset are named by the stem of its input file, which for a retrieved dataset is
`<accession>_<Author>_<year>`: `R_T_vs_A_H_segmented_<accession>_<Author>_<year>.csv`. The label
stays `Author year`, and [`curve_accessions`](@ref) gives the accession of each.
"""
function write_results(result::ExtractionResult, directory::AbstractString)
    isdir(directory) &&
        !isempty(readdir(directory)) &&
        throw(
            ArgumentError(
                "$(directory) already holds files; write_results never writes into a directory \
                 that does, move it aside first"
            ),
        )
    configuration = result.configuration
    identifier = run_identifier(configuration)
    # The files named by the identifier; a name longer than a file system admits would fail
    # halfway through writing the run, so it is refused before anything is written.
    for name in (
        "manifest_$(identifier).toml",
        "segmented_curves_$(identifier).csv",
        "total_average_R_T_$(identifier).csv",
        "leave_one_out_$(identifier).csv",
    )
        ncodeunits(name) <= MAX_FILE_NAME_BYTES || throw(
            ArgumentError("the file name $(name) is $(ncodeunits(name)) bytes, beyond the \
                 $(MAX_FILE_NAME_BYTES) a file system admits")
        )
    end
    mkpath(directory)
    digits = configuration.output.significant_digits
    written = Dict{String,String}()
    tokens = _file_tokens(result)
    accessions = curve_accessions(result)

    for (curves, name, quantity) in (
        (result.r_ν, "r_nu_vs_A_H", "r_nu"),
        (result.R_T, "R_T_vs_A_H", "R_T"),
        ([result.consensus_r_ν], "r_nu_vs_A_H_consensus", "r_nu"),
    )
        for curve in curves
            isempty(curve) && continue
            path = joinpath(directory, "$(name)_$(tokens[curve.label]).csv")
            CSV.write(path, _ratio_table(curve, quantity, digits))
            written["$(name)/$(curve.label)"] = path
        end
    end

    for curve in result.segmented_curves
        token = tokens[curve.label]
        points = pivots(curve.fit)
        pivot_table = DataFrame(;
            A_H = [point[1] for point in points],
            r_nu = [round(point[2]; sigdigits = digits) for point in points],
            # The propagated uncertainty at each pivot. The fit carries it, and a file that states
            # the parameterization without it invites the reader to assume there is none.
            r_nu_uncertainty = [
                round(last(evaluate(curve.fit, Float64(point[1]))); sigdigits = digits) for
                point in points
            ],
        )
        path = joinpath(directory, "r_nu_vs_A_H_pivots_$(token).csv")
        CSV.write(path, pivot_table)
        written["r_nu_pivots/$(curve.label)"] = path

        for (segmented, name, quantity) in (
            (curve.r_ν, "r_nu_vs_A_H_segmented", "r_nu"),
            (curve.R_T, "R_T_vs_A_H_segmented", "R_T"),
        )
            path = joinpath(directory, "$(name)_$(token).csv")
            CSV.write(path, _ratio_table(segmented, quantity, digits))
            written["$(quantity)_segmented/$(curve.label)"] = path
        end
    end

    # What a run knows about each segmented curve beyond the manifest: one row per manifest curve,
    # keyed by its label. No column is an energy; the ratios are dimensionless.
    excluded = _excluded_labels(
        configuration.excluded_datasets, result.datasets, "multiplicity.exclude"
    )
    pooled = Set(pooled_datasets(result))
    # The other datasets of its experiment, for a dataset pooled as one with them.
    partners = Dict(
        label => filter(!=(label), group) for group in result.correlation_groups for
        label in group
    )
    autocorrelation = deviation_autocorrelation(result)
    rows = [
        (
            label = curve.label,
            accession = _text_field(accessions[curve.label]),
            kind = curve.kind,
            pooled = curve.kind == "systematic_trend" || curve.label in pooled,
            segments = segments(curve.fit),
            pinned_at_symmetric_split = curve.fit.pinned_value !== nothing,
            first_A_H = first(curve.R_T.A_H),
            last_A_H = last(curve.R_T.A_H),
            pairs = curve.pairs,
            measured_points = round(curve.fit.measured_points; sigdigits = digits),
            coverage = round(curve.coverage; sigdigits = digits),
            reduced_chi_squared = round(curve.fit.wrss / curve.fit.dof; sigdigits = digits),
            chi_squared_over_expectation = round(
                curve.fit.wrss / curve.fit.expected_wrss; sigdigits = digits
            ),
            weights_imputed = curve.fit.weights_imputed,
            range_mean_R_T = round(
                first(result.range_mean_R_T[curve.label]); sigdigits = digits
            ),
            range_mean_R_T_uncertainty = round(
                last(result.range_mean_R_T[curve.label]); sigdigits = digits
            ),
            # The autocorrelation the covariance of the trend was formed with.
            deviation_autocorrelation = if curve.kind == "systematic_trend"
                if ismissing(autocorrelation)
                    missing
                else
                    round(autocorrelation; sigdigits = digits)
                end
            else
                missing
            end,
        ) for curve in result.segmented_curves
    ]
    path = joinpath(directory, "segmented_curves_$(identifier).csv")
    CSV.write(path, DataFrame(rows))
    written["segmented_curves"] = path

    # The total average over each yield distribution: one row per segmented curve and
    # distribution, which is how the literature tabulates it. The covariance-propagated
    # uncertainty first; the independent-points one, the published approximation, beside it; then
    # the fraction of the distribution's yield the curve takes in. For the systematic trend,
    # ⟨R_T⟩ refitted with one and two segments more than selected, and the jackknife uncertainty
    # over the pooled datasets with the extremes of ⟨R_T⟩ with one of them left out.
    if !isempty(result.total_average_R_T)
        rows = NamedTuple{
            (
                :segmented_curve,
                :mass_yield,
                :R_T,
                :R_T_uncertainty,
                :R_T_uncertainty_independent_points,
                :R_T_uncertainty_leave_one_out,
                :R_T_leave_one_out_min,
                :R_T_leave_one_out_max,
                :yield_fraction,
                :R_T_one_more_segment,
                :R_T_two_more_segments,
            ),
            Tuple{
                String,
                String,
                Float64,
                Float64,
                Float64,
                Union{Missing,Float64},
                Union{Missing,Float64},
                Union{Missing,Float64},
                Float64,
                Union{Missing,Float64},
                Union{Missing,Float64},
            },
        }[]
        for curve in result.segmented_curves
            per_distribution = get(result.total_average_R_T, curve.label, nothing)
            per_distribution === nothing && continue
            for distribution in result.mass_yields
                entry = get(per_distribution, distribution.label, nothing)
                entry === nothing && continue
                more = if curve.kind == "systematic_trend"
                    get(result.segment_count_sensitivity, distribution.label, Tuple{Int,Float64}[])
                else
                    Tuple{Int,Float64}[]
                end
                alternative(extra) = begin
                    found = findfirst(t -> t[1] == segments(curve.fit) + extra, more)
                    found === nothing ? missing : round(more[found][2]; sigdigits = digits)
                end
                spread = if curve.kind == "systematic_trend"
                    leave_one_out_spread(result, distribution.label)
                else
                    nothing
                end
                spread_field(field) =
                    spread === nothing ? missing : round(spread[field]; sigdigits = digits)
                push!(
                    rows,
                    (
                        segmented_curve = curve.label,
                        mass_yield = distribution.label,
                        R_T = round(entry.value; sigdigits = digits),
                        R_T_uncertainty = round(entry.uncertainty; sigdigits = digits),
                        R_T_uncertainty_independent_points = round(
                            entry.uncertainty_independent_points; sigdigits = digits
                        ),
                        R_T_uncertainty_leave_one_out = spread_field(:uncertainty),
                        R_T_leave_one_out_min = spread_field(:min),
                        R_T_leave_one_out_max = spread_field(:max),
                        yield_fraction = round(entry.yield_fraction; sigdigits = digits),
                        R_T_one_more_segment = alternative(1),
                        R_T_two_more_segments = alternative(2),
                    ),
                )
            end
        end
        if !isempty(rows)
            path = joinpath(directory, "total_average_R_T_$(identifier).csv")
            CSV.write(path, DataFrame(rows))
            written["total_average_R_T"] = path
        end
    end

    # The systematic trend refitted with each pooled dataset left out: one row per dataset and
    # yield distribution it has ⟨R_T⟩ over, a single row without a distribution where the refit
    # has none, so that every refit, a failed one included, is on record.
    if !isempty(result.leave_one_out)
        rows = NamedTuple{
            (
                :dataset_left_out,
                :accession,
                :mass_yield,
                :R_T,
                :segments,
                :breakpoints,
                :reduced_chi_squared,
                :outcome,
            ),
            Tuple{
                String,
                Union{Missing,String},
                Union{Missing,String},
                Union{Missing,Float64},
                Union{Missing,Int},
                Union{Missing,String},
                Union{Missing,Float64},
                String,
            },
        }[]
        for entry in result.leave_one_out
            common = (
                dataset_left_out = join(entry.datasets, " + "),
                accession = _text_field(
                    join(
                        filter(!isempty, [get(accessions, l, "") for l in entry.datasets]), " "
                    ),
                ),
            )
            refit = (
                segments = entry.segments == 0 ? missing : entry.segments,
                breakpoints = _text_field(join(entry.breakpoints, " ")),
                reduced_chi_squared = if isnan(entry.reduced_chi_squared)
                    missing
                else
                    round(entry.reduced_chi_squared; sigdigits = digits)
                end,
                outcome = entry.outcome,
            )
            if isempty(entry.total_average_R_T)
                push!(rows, (; common..., mass_yield = missing, R_T = missing, refit...))
                continue
            end
            for distribution in result.mass_yields
                value = get(entry.total_average_R_T, distribution.label, nothing)
                value === nothing && continue
                push!(
                    rows,
                    (;
                        common...,
                        mass_yield = distribution.label,
                        R_T = round(value; sigdigits = digits),
                        refit...,
                    ),
                )
            end
        end
        path = joinpath(directory, "leave_one_out_$(identifier).csv")
        CSV.write(path, DataFrame(rows))
        written["leave_one_out"] = path
    end

    # Per-dataset diagnostics, written for every dataset whether or not it was pooled or fitted,
    # with the reaction-code qualifiers the retrieval recorded for it.
    if !isempty(result.dataset_diagnostics)
        rows = [
            (
                dataset = d.label,
                accession = _text_field(accessions[d.label]),
                points = d.points,
                pairs = d.pairs,
                first_pair = d.first_pair,
                last_pair = d.last_pair,
                coverage = round(d.coverage; sigdigits = digits),
                outside_physical_range = d.outside_physical_range,
                symmetry_departure = d.symmetry_departure,
                complement_sum = d.complement_sum,
                complement_spread = d.complement_spread,
                without_uncertainties = d.without_uncertainties,
                qualifiers = _text_field(join(get(result.qualifiers, d.label, String[]), " ")),
                pooling_weight = round(
                    pooling_weight(retrieval_record(result.datasets[i].source));
                    sigdigits = digits,
                ),
                flagged = !isempty(_flagged(get(result.qualifiers, d.label, String[]))),
                pair_sum_deviation = _scale_field(result.datasets[i], :deviation, digits),
                scale_consistent = _scale_field(result.datasets[i], :consistent, digits),
                pooled = d.label in pooled,
                pooled_with = _text_field(join(get(partners, d.label, String[]), " + ")),
                exclusion_reason = _text_field(get(excluded, d.label, "")),
                segmented_curve = _text_field(get(result.dataset_outcomes, d.label, "")),
                curve_flagged = haskey(result.curve_flags, d.label),
                curve_flag_reason = _text_field(get(result.curve_flags, d.label, "")),
            ) for (i, d) in enumerate(result.dataset_diagnostics)
        ]
        path = joinpath(directory, "dataset_diagnostics.csv")
        CSV.write(path, DataFrame(rows))
        written["dataset_diagnostics"] = path
    end

    written["manifest"] = _write_manifest(result, directory, identifier, written)
    @info "results written" directory
    return written
end

# The longest file name, in bytes, that common file systems admit (NAME_MAX).
const MAX_FILE_NAME_BYTES = 255

_file_token(label::AbstractString) = replace(strip(label), r"[^A-Za-z0-9.\-]+" => "_")

# A text field of a table: one that holds nothing is written as an empty field, as an absent
# number is, and not as a quoted empty string.
_text_field(text::AbstractString) = isempty(text) ? missing : String(text)
