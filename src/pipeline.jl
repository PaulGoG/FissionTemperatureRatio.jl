# The pipeline: from experimental multiplicity data to a segmented temperature ratio.

"""
    SymmetryDiagnostics

The identities that hold by construction at the symmetric split, checked rather than assumed.

# Fields

- `charge_set_invariant`: whether the charge numbers retained at the symmetric split are
  invariant under `Z -> Z₀ - Z`, which is what makes `R_a = 1` exact there; `missing` for a
  fissioning nucleus of odd mass number, which has no symmetric split.
- `R_a_at_symmetric_split`: the level density parameter ratio there, one by identity; `missing`
  where there is no symmetric split or the ratio is undefined.
- `pinned_curves`: the labels of the segmented curves pinned to `r_ν = 1/2` at the symmetric
  split.
- `warnings`: every departure found, one sentence each; empty when the identities hold.
"""
struct SymmetryDiagnostics
    charge_set_invariant::Union{Bool,Missing}
    R_a_at_symmetric_split::Union{Float64,Missing}
    pinned_curves::Vector{String}
    warnings::Vector{String}
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
- `R_a`: the level density parameter ratio against heavy-fragment mass number.
- `segmented_curves`: one per dataset that supports a fit and reaches the coverage floor,
  followed by the systematic-trend curve. Alternatives offered to a prompt emission code, not an
  ensemble.
- `range_mean_R_T`: for each segmented curve, the mean of its temperature ratio over its mass
  range with the uncertainty propagated through the curve's covariance; see [`range_mean`](@ref)
  for why this is not the total average.
- `mass_yields`: the fragment mass yield distributions read, empty when the configuration names
  none.
- `total_average_R_T`: the quantity the literature quotes, `⟨R_T⟩` over each yield distribution,
  keyed by curve label and then by yield label. Empty when no yields were given.
- `consensus_r_ν`: the combined multiplicity ratio the systematic-trend curve was fitted to, so
  that the trend can be checked against its own input rather than taken on trust.
- `dataset_diagnostics`: one record per multiplicity dataset read, in the order read.
- `dataset_outcomes`: for every dataset, whether it offers a segmented curve and, if not, why.
- `symmetry`: the identities at the symmetric split, checked.
"""
struct ExtractionResult
    configuration::Configuration
    datasets::Vector{Multiplicity}
    r_ν::Vector{RatioCurve}
    R_T::Vector{RatioCurve}
    R_a::Dict{Int,Float64}
    segmented_curves::Vector{SegmentedCurve}
    range_mean_R_T::Dict{String,Tuple{Float64,Float64}}
    mass_yields::Vector{MassYield}
    total_average_R_T::Dict{String,Dict{String,TotalAverage}}
    consensus_r_ν::RatioCurve
    dataset_diagnostics::Vector{DatasetDiagnostics}
    dataset_outcomes::Dict{String,String}
    symmetry::SymmetryDiagnostics
end

"""
    systematic_trend(result) -> SegmentedCurve

The segmented curve that follows the systematic behaviour of the ratio rather than any single
dataset.
"""
function systematic_trend(result::ExtractionResult)
    index = findfirst(c -> c.label == TREND_LABEL, result.segmented_curves)
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
"""
function read_multiplicity_directory(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("multiplicity directory not found: $(directory)"))
    files = _data_files(directory)
    isempty(files) &&
        throw(ArgumentError("multiplicity directory holds no data files: $(directory)"))
    labels = _dataset_label.(files)
    for (index, file) in enumerate(files)
        count(==(labels[index]), labels) == 1 && continue
        labels[index] = "$(labels[index]) ($(_accession(file)))"
    end
    return [
        read_multiplicity(joinpath(directory, file); label = label) for
        (file, label) in zip(files, labels)
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
    run_pipeline(configuration) -> ExtractionResult

Execute the extraction for one configuration.

The steps are: read the tabulated data and the experimental multiplicity data; build the
fragmentation range and the isobaric charge distribution; form the multiplicity ratio
`r_ν = ν_H/(ν_L + ν_H)` of every dataset and the temperature ratio from it; and describe the
multiplicity ratio by joined straight segments, once per dataset and once more following its
systematic behaviour.

A run therefore returns several segmented curves rather than one. Where the datasets of a
fissioning nucleus disagree — as they do markedly for some — a curve fitted through all of them
describes none of them, and a prompt emission code has nothing to validate it against. They are
offered as alternatives: the code takes one as input, and which one describes reality is settled
by comparing the multiplicity distributions and yields it produces against experiment. The
systematic-trend curve is the one to use where a dataset is too sparse or too scattered to
determine the shape on its own.

The temperature ratio is obtained by transforming the fitted multiplicity ratio, not by fitting
segments to the temperature ratio a second time. The transformation is exact, whereas a second fit
would discard the uncertainty of the first and impose a piecewise-linear shape on a quantity that
is not piecewise-linear.

Nothing is written. [`write_results`](@ref) writes the tables and the manifest of a result into a
directory of the caller's choosing; `scripts/run.jl` names that directory by the run identifier
and adds the provenance record and the figures.
"""
function run_pipeline(configuration::Configuration)
    @info "reading input" configuration = configuration.source
    model = build_level_density_model(configuration.level_density)
    charge_data = read_charge_distribution(
        configuration.fragmentation.charge_distribution_file,
        configuration.fragmentation.fallback_charge_polarization,
        configuration.fragmentation.fallback_charge_dispersion,
    )
    datasets = read_multiplicity_directory(configuration.multiplicity_directory)

    A₀ = configuration.system.A₀
    Z₀ = configuration.system.Z₀
    range = A_H_range(configuration)

    domain = fragmentation_domain(
        A₀, Z₀, range, configuration.fragmentation.charges_per_mass, charge_data
    )
    if !isempty(domain.fallback_masses)
        @warn "tabulated charge polarization and dispersion unavailable; the average values were \
               used" masses = _summarize_range(domain.fallback_masses) ΔZ =
            configuration.fragmentation.fallback_charge_polarization σ_Z =
            configuration.fragmentation.fallback_charge_dispersion
    end

    R_a = level_density_ratio(
        configuration.level_density.ratio_averaging, model, A₀, Z₀, domain
    )
    isempty(R_a) &&
        throw(ArgumentError("the level density parameter ratio is undefined over the whole \
                       fragmentation range; check the mass excess table"))

    r_ν = [multiplicity_ratio(data, A₀, range) for data in datasets]
    R_T = [temperature_ratio(curve, R_a) for curve in r_ν]
    diagnostics = [diagnose(datasets[i], r_ν[i], A₀, range) for i in eachindex(datasets)]
    usable = findall(!isempty, r_ν)
    isempty(usable) &&
        throw(ArgumentError("no dataset provides both fragments of any pair within \
                       $(first(range)):$(last(range))"))

    settings = configuration.segments
    segmented_curves = SegmentedCurve[]
    outcomes = Dict{String,String}()

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
        if coverage < settings.min_dataset_coverage
            outcomes[data.label] = "coverage $(round(coverage; digits = 3)) below the floor \
                                    $(settings.min_dataset_coverage); pooled only"
            @info "no segmented curve for this dataset" dataset = data.label coverage floor =
                settings.min_dataset_coverage
            continue
        end
        windows =
            settings.windows_apply_to_datasets ? settings.required_windows : UnitRange{Int}[]
        segmented = _segment(curve, R_a, settings, curve.label, windows, A₀)
        if segmented isa SegmentedCurve
            push!(segmented_curves, segmented)
            outcomes[data.label] = "segmented curve"
        else
            outcomes[data.label] = "no fit: $(segmented)"
        end
    end
    unpinned = [c.label for c in segmented_curves if c.fit.pinned_value === nothing]
    if settings.pin_symmetric_split && !isempty(unpinned)
        @info "fitted without the pin: no complete pair at the symmetric split" datasets =
            unpinned
    end

    # One curve following the systematic behaviour of the ratio: the minimum at the heavy magic
    # fragment placed rather than fitted, and the rise above the most probable fragmentation taken
    # through the whole body of data, so its slope falls between those of the individual datasets.
    # This is the curve to use where a dataset is too sparse or too scattered to determine the
    # shape on its own.
    # Datasets named in the configuration are kept out of the pooling but not out of the run: they
    # are still fitted, written and diagnosed, so an exclusion is visible rather than a silent
    # absence.
    admitted = [
        i for i in usable if !haskey(configuration.excluded_datasets, datasets[i].label)
    ]
    for (label, reason) in configuration.excluded_datasets
        any(data.label == label for data in datasets) ||
            @warn "configuration excludes a dataset that was not read" dataset = label reason
    end
    isempty(admitted) && throw(
        ArgumentError("every usable dataset is excluded from the pooling by configuration")
    )
    pooled = consensus(r_ν[admitted]; label = TREND_LABEL)
    trend = _segment(pooled, R_a, settings, TREND_LABEL, settings.required_windows, A₀)
    if !(trend isa SegmentedCurve) && !isempty(settings.required_windows)
        @warn "no segmented curve satisfies the required windows; the systematic-trend curve \
               was fitted without them" windows = settings.required_windows reason = trend
        trend = _segment(pooled, R_a, settings, TREND_LABEL, UnitRange{Int}[], A₀)
    end
    trend isa SegmentedCurve && push!(segmented_curves, trend)

    isempty(segmented_curves) &&
        throw(ArgumentError("no dataset supports a segmented curve with \
                       min_points_per_segment = $(settings.min_points_per_segment)"))

    for curve in segmented_curves
        @info "segmented curve" dataset = curve.label segments = segments(curve.fit) breakpoints =
            curve.fit.breakpoints reduced_chi_squared = curve.fit.wrss / curve.fit.dof
    end

    symmetry = _symmetry_diagnostics(configuration, domain, R_a, segmented_curves)
    for message in symmetry.warnings
        @warn message
    end

    mass_yields = if configuration.yield_directory === nothing
        MassYield[]
    else
        read_mass_yield_directory(configuration.yield_directory)
    end
    total_average_R_T = _total_averages(segmented_curves, mass_yields)
    for curve in segmented_curves
        for distribution in mass_yields
            entry = get(
                get(total_average_R_T, curve.label, Dict()), distribution.label, nothing
            )
            entry === nothing && continue
            @info "total average" dataset = curve.label yield = distribution.label R_T =
                entry.value uncertainty = entry.uncertainty
        end
    end

    return ExtractionResult(
        configuration,
        datasets,
        r_ν,
        R_T,
        R_a,
        segmented_curves,
        Dict(c.label => range_mean(c) for c in segmented_curves),
        mass_yields,
        total_average_R_T,
        pooled,
        diagnostics,
        outcomes,
        symmetry,
    )
end

# Whether a yield distribution carries positive weight at any mass number of the curve.
function _overlaps(curve::RatioCurve, distribution::MassYield)
    total = 0.0
    for mass in curve.A_H
        entry = mass_yield(distribution, mass)
        ismissing(entry) || (total += entry[1])
    end
    return total > 0
end

# The total average of every segmented curve over every yield distribution. A pair that shares no
# mass number is omitted rather than reported as zero: a distribution covering only the light wing
# says nothing about a ratio defined on the heavy one.
function _total_averages(curves::Vector{SegmentedCurve}, mass_yields::Vector{MassYield})
    averages = Dict{String,Dict{String,TotalAverage}}()
    isempty(mass_yields) && return averages
    for curve in curves
        per_distribution = Dict{String,TotalAverage}()
        for distribution in mass_yields
            if _overlaps(curve.R_T, distribution)
                per_distribution[distribution.label] = TotalAverage(curve, distribution)
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
# few mass pairs.
#
# `fit_segments` pins at the first abscissa it is given. r_ν = 1/2 is an identity at A₀/2 only, so
# the pin is applied to a curve whose first complete pair is the symmetric split and to no other;
# a dataset starting above it is fitted unpinned.
function _segment(
    curve::RatioCurve,
    R_a::AbstractDict{Int,Float64},
    settings::SegmentSettings,
    label::AbstractString,
    windows::Vector{UnitRange{Int}},
    A₀::Integer,
)
    pinned = settings.pin_symmetric_split && !isempty(curve) && 2 * first(curve.A_H) == A₀
    fit = try
        fit_segments(
            curve.A_H,
            curve.ratio,
            curve.σ;
            max_segments = settings.max_segments,
            min_points_per_segment = settings.min_points_per_segment,
            min_segment_span = settings.min_segment_span,
            pinned_value = pinned ? 0.5 : nothing,
            required_windows = windows,
            bounds = (0.0, 1.0),
        )
    catch exception
        exception isa InsufficientDataError || rethrow()
        @warn "no segmented curve for this dataset" dataset = label reason = exception.msg
        return exception.msg
    end
    return SegmentedCurve(label, fit, R_a)
end

# The identities that hold by construction at the symmetric split, checked rather than assumed.
function _symmetry_diagnostics(
    configuration::Configuration,
    domain::FragmentationDomain,
    R_a::AbstractDict{Int,Float64},
    curves::Vector{SegmentedCurve},
)
    warnings = String[]
    A₀ = configuration.system.A₀
    invariant = symmetric_charge_set_is_invariant(domain, A₀, configuration.system.Z₀)
    if invariant === false
        push!(
            warnings,
            "the charge numbers retained at the symmetric split are not invariant under \
             Z -> Z₀ - Z, so R_a = 1 and R_T = 1 hold there only approximately; this occurs when \
             Z₀ is odd and the most probable charge falls between two integers",
        )
    end

    R_a_symmetric = missing
    if iseven(A₀) && haskey(R_a, A₀ ÷ 2)
        R_a_symmetric = R_a[A₀ ÷ 2]
        deviation = abs(R_a_symmetric - 1)
        deviation < 1e-8 || push!(
            warnings,
            "R_a departs from unity by $(round(deviation; sigdigits = 3)) at the symmetric \
             split, where the two fragments are the same nuclide and the ratio must be exactly one",
        )
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
    return SymmetryDiagnostics(invariant, R_a_symmetric, pinned, warnings)
end

# A machine-readable index of what a run produced, for a code that consumes these curves rather
# than a person reading them. It names the system, lists every segmented curve with the file to
# read for it, and says which datasets were pooled. Deliberately separate from the run metadata,
# which records how the result was produced; this records what is on offer.
function _write_manifest(
    result::ExtractionResult,
    directory::AbstractString,
    identifier::AbstractString,
    written::Dict{String,String},
)
    entries = Dict{String,Any}[]
    coverage = Dict(d.label => d.coverage for d in result.dataset_diagnostics)
    range = A_H_range(result.configuration)
    for curve in result.segmented_curves
        label = curve.label
        trend = label == TREND_LABEL
        # The pairs the curve was fitted to: the dataset's own, or the combined curve's.
        fitted_to = if trend
            result.consensus_r_ν
        else
            result.r_ν[findfirst(c -> c.label == label, result.r_ν)]
        end
        averages = get(result.total_average_R_T, label, Dict{String,TotalAverage}())
        entry = Dict{String,Any}(
            "label" => label,
            "kind" => trend ? "systematic_trend" : "dataset",
            "pooled" => trend || !haskey(result.configuration.excluded_datasets, label),
            "segments" => segments(curve.fit),
            "pinned_at_symmetric_split" => curve.fit.pinned_value !== nothing,
            "breakpoints" => curve.fit.breakpoints,
            # Parallel arrays rather than pairs: a mixed array of integers and floats is
            # promoted to floats on serialization, and a mass number is not a float.
            "pivot_A_H" => [p[1] for p in pivots(curve.fit)],
            "pivot_r_nu" => [p[2] for p in pivots(curve.fit)],
            "reduced_chi_squared" => curve.fit.wrss / curve.fit.dof,
            # What the curve rests on: the pairs it was fitted to and the span they cover, so a
            # consumer can see a sparse curve for what it is without reading the data.
            "first_A_H" => first(curve.R_T.A_H),
            "last_A_H" => last(curve.R_T.A_H),
            "pairs" => length(fitted_to),
            "coverage" =>
                trend ? count(in(range), fitted_to.A_H) / length(range) : coverage[label],
            "range_mean_R_T" => collect(result.range_mean_R_T[label]),
            "total_average_R_T" => Dict{String,Any}(
                k => [v.value, v.uncertainty, v.uncertainty_independent_points] for
                (k, v) in averages
            ),
        )
        # The file a consumer should read. The temperature ratio is tabulated at every mass
        # number, because it is not piecewise-linear even where the multiplicity ratio is:
        # interpolating between the segment pivots would cut across the structure the level
        # density parameter ratio puts into it.
        for (key, name) in (
            ("R_T_segmented/$(label)", "temperature_ratio_file"),
            ("r_nu_pivots/$(label)", "multiplicity_ratio_pivots_file"),
        )
            haskey(written, key) && (entry[name] = basename(written[key]))
        end
        push!(entries, entry)
    end

    manifest = Dict{String,Any}(
        "system" => _system_record(result.configuration.system),
        "run" => Dict{String,Any}(
            "identifier" => identifier,
            "package_version" => string(PACKAGE_VERSION),
            "quantity" => "R_T = T_L/T_H of complementary fully accelerated fragments",
            "ordinate" => "R_T",
            "abscissa" => ["A_H"],
            # Stated so that a reader knows what the columns hold, not so that a reader looks
            # them up by name: the files are read by column position.
            "columns" => ["A_H", "R_T", "R_T_uncertainty"],
            "total_average_R_T_columns" =>
                ["R_T", "R_T_uncertainty", "R_T_uncertainty_independent_points"],
        ),
        "segmented_curve" => entries,
    )
    path = joinpath(directory, "manifest_$(identifier).toml")
    open(path, "w") do io
        return TOML.print(io, manifest; sorted = true)
    end
    return path
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

Write the tabulated ratios, the segment pivots, the total averages, the dataset diagnostics and
the manifest of a run into `directory`, and return the paths written, keyed by content.

`directory` is created. One that already holds files is refused rather than written into, so a
caller that wants a second run of one configuration beside the first moves the first aside;
`scripts/run.jl` does, numbering it `#1`, `#2`, …, and names the directory
`data/sims/<system>/<run identifier>/`. The provenance record and the figures are not written
here: [`run_metadata`](@ref) supplies what the library knows about a run, and the script adds the
rest, writes `metadata.toml`, and calls [`write_figures`](@ref).

The manifest is `manifest_<run identifier>.toml`, the one file whose name carries the identifier:
a consuming code stages the whole directory and selects the manifest by that prefix.
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
    mkpath(directory)
    configuration = result.configuration
    identifier = run_identifier(configuration)
    digits = configuration.output.significant_digits
    written = Dict{String,String}()

    for (curves, name, quantity) in (
        (result.r_ν, "r_nu_vs_A_H", "r_nu"),
        (result.R_T, "R_T_vs_A_H", "R_T"),
        ([result.consensus_r_ν], "r_nu_vs_A_H_consensus", "r_nu"),
    )
        for curve in curves
            isempty(curve) && continue
            path = joinpath(directory, "$(name)_$(_file_token(curve.label)).csv")
            CSV.write(path, _ratio_table(curve, quantity, digits))
            written["$(name)/$(curve.label)"] = path
        end
    end

    for curve in result.segmented_curves
        token = _file_token(curve.label)
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

    # The total average over each yield distribution: one row per segmented curve and
    # distribution, which is how the literature tabulates it. The covariance-propagated
    # uncertainty first; the independent-points one, the published approximation, beside it.
    if !isempty(result.total_average_R_T)
        rows = NamedTuple{
            (
                :segmented_curve,
                :mass_yield,
                :R_T,
                :R_T_uncertainty,
                :R_T_uncertainty_independent_points,
            ),
            Tuple{String,String,Float64,Float64,Float64},
        }[]
        for curve in result.segmented_curves
            per_distribution = get(result.total_average_R_T, curve.label, nothing)
            per_distribution === nothing && continue
            for distribution in result.mass_yields
                entry = get(per_distribution, distribution.label, nothing)
                entry === nothing && continue
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
                    ),
                )
            end
        end
        if !isempty(rows)
            path = joinpath(directory, "total_average_R_T.csv")
            CSV.write(path, DataFrame(rows))
            written["total_average_R_T"] = path
        end
    end

    # Per-dataset diagnostics, written for every dataset whether or not it was pooled or fitted.
    if !isempty(result.dataset_diagnostics)
        excluded = configuration.excluded_datasets
        rows = [
            (
                dataset = d.label,
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
                pooled = !haskey(excluded, d.label),
                exclusion_reason = get(excluded, d.label, ""),
                segmented_curve = get(result.dataset_outcomes, d.label, ""),
            ) for d in result.dataset_diagnostics
        ]
        path = joinpath(directory, "dataset_diagnostics.csv")
        CSV.write(path, DataFrame(rows))
        written["dataset_diagnostics"] = path
    end

    written["manifest"] = _write_manifest(result, directory, identifier, written)
    @info "results written" directory
    return written
end

_file_token(label::AbstractString) = replace(strip(label), r"[^A-Za-z0-9.\-]+" => "_")
