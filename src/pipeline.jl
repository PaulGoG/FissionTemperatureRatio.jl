# The pipeline: from experimental multiplicity data to a segmented temperature ratio.

"""
    SegmentedCurve

One piecewise-linear description of the multiplicity ratio, and the temperature ratio obtained
from it.

A run produces several: one per experimental dataset, and one following the systematic behaviour
of the ratio. They are alternatives, not an ensemble to be averaged — a prompt emission code takes
one of them as input, and which one describes reality is settled downstream, by comparing the
multiplicity distributions and yields that code produces against experiment.

# Fields

- `label`: the dataset the curve came from, or `"systematic trend"`.
- `fit`: the piecewise-linear fit of `r_ν`.
- `r_ν`, `R_T`: the segmented multiplicity ratio and the temperature ratio from it.
"""
struct SegmentedCurve
    label::String
    fit::SegmentedFit
    r_ν::RatioCurve
    R_T::RatioCurve
end

"""
    ExtractionResult

Everything a run produces, held together so that it can be inspected interactively as well as
written to disk.

# Fields

- `configuration`: the configuration the run was driven by.
- `datasets`: the experimental multiplicity data that was read.
- `r_ν`, `R_T`: the multiplicity and temperature ratios extracted from each dataset.
- `R_a`: the level density parameter ratio against heavy-fragment mass number.
- `segmented_curves`: one per dataset that supports a fit, followed by the systematic-trend curve.
  Alternatives offered to a prompt emission code, not an ensemble.
- `range_mean_R_T`: for each segmented curve, the inverse-variance weighted mean of its
  temperature ratio over the fragment mass range, with its uncertainty. This is *not* the total
  average quoted in the literature: it weights every mass number equally and so is dominated by
  the far-asymmetric tail, where the yield is negligible.
- `mass_yields`: the fragment mass yield distributions read, empty when the configuration names
  none.
- `consensus_r_ν`: the combined multiplicity ratio the systematic-trend curve was fitted to, so
  that the trend can be checked against its own input rather than taken on trust.
- `dataset_diagnostics`: one record per multiplicity dataset read, in the order read.
- `total_average_R_T`: the quantity the literature quotes, `⟨R_T⟩` over each yield distribution,
  keyed by curve label and then by yield label. Empty when no yields were given.
- `diagnostics`: checks of the exact identities at the symmetric split.
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
    total_average_R_T::Dict{String,Dict{String,Tuple{Float64,Float64}}}
    consensus_r_ν::RatioCurve
    dataset_diagnostics::Vector{DatasetDiagnostics}
    diagnostics::Dict{String,Any}
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
literature.
"""
function read_multiplicity_directory(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("multiplicity directory not found: $(directory)"))
    files = _data_files(directory)
    isempty(files) &&
        throw(ArgumentError("multiplicity directory holds no data files: $(directory)"))
    return [
        read_multiplicity(joinpath(directory, file); label = _dataset_label(file)) for
        file in files
    ]
end

# Data files carry the `.dat` extension; anything else in the directory — a run record, a note —
# is not data and is skipped rather than parsed and rejected.
function _data_files(directory::AbstractString)
    return sort!(filter!(f -> endswith(f, ".dat") && !startswith(f, "."), readdir(directory)))
end

# `<accession>_<Author>_<year>.dat`.
function _dataset_label(file::AbstractString)
    stem = replace(splitext(file)[1], r"^[0-9]+_" => "")
    label = replace(stem, '_' => ' ')
    # Initials are written without a space in the file names; restore it for the legend.
    return replace(label, r"(?<=\b[A-Z])\.(?=[A-Z][a-z])" => ". ")
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
    σ = reduce(vcat, (curve.σ for curve in curves); init = Float64[])
    order = sortperm(A_H)
    return RatioCurve(A_H[order], ratio[order], σ[order], String(label))
end

"""
    run_pipeline(configuration; write_output = true) -> ExtractionResult

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

With `write_output`, tabulated results, figures and run metadata are written under `output_root`;
an existing file is never overwritten, so a rerun cannot destroy a previous result.
"""
function run_pipeline(
    configuration::Configuration;
    write_output::Bool = true,
    output_root::AbstractString = projectdir(),
)
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
    usable = findall(!isempty, r_ν)
    isempty(usable) &&
        throw(ArgumentError("no dataset provides both fragments of any pair within \
                       $(first(range)):$(last(range))"))

    settings = configuration.segments
    segmented_curves = SegmentedCurve[]

    # One curve per dataset. These are the alternatives a prompt emission code chooses between.
    for index in usable
        curve = r_ν[index]
        windows =
            settings.windows_apply_to_datasets ? settings.required_windows : UnitRange{Int}[]
        segmented = _segment(curve, R_a, settings, curve.label, windows, A₀)
        segmented === nothing && continue
        push!(segmented_curves, segmented)
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
            @warn "configuration excludes a dataset that was not read" dataset = label
    end
    isempty(admitted) && throw(
        ArgumentError("every usable dataset is excluded from the pooling by configuration")
    )
    pooled = consensus(r_ν[admitted]; label = TREND_LABEL)
    trend = _segment(pooled, R_a, settings, TREND_LABEL, settings.required_windows, A₀)
    if trend === nothing && !isempty(settings.required_windows)
        @warn "no segmented curve satisfies the required windows; the systematic-trend curve \
               was fitted without them" windows = settings.required_windows
        trend = _segment(pooled, R_a, settings, TREND_LABEL, UnitRange{Int}[], A₀)
    end
    trend === nothing || push!(segmented_curves, trend)

    isempty(segmented_curves) &&
        throw(ArgumentError("no dataset supports a segmented curve with \
                       min_points_per_segment = $(settings.min_points_per_segment)"))

    for curve in segmented_curves
        @info "segmented curve" dataset = curve.label segments = segments(curve.fit) breakpoints =
            curve.fit.breakpoints reduced_chi_squared = curve.fit.wrss / curve.fit.dof
    end

    diagnostics = _symmetry_diagnostics(configuration, domain, R_a, segmented_curves)
    for (key, message) in diagnostics["warnings"]
        @warn message identity = key
    end

    diagnostics_by_dataset = [diagnose(datasets[i], r_ν[i], A₀) for i in eachindex(datasets)]

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
            @info "total average" dataset = curve.label yield = distribution.label R_T = entry[1] uncertainty = entry[2]
        end
    end

    result = ExtractionResult(
        configuration,
        datasets,
        r_ν,
        R_T,
        R_a,
        segmented_curves,
        Dict(c.label => weighted_mean(c.R_T) for c in segmented_curves),
        mass_yields,
        total_average_R_T,
        pooled,
        diagnostics_by_dataset,
        diagnostics,
    )

    write_output && write_results(result; root = output_root)
    return result
end

# The total average of every segmented curve over every yield distribution. A pair that shares no
# mass number is omitted rather than reported as zero: a distribution covering only the light wing
# says nothing about a ratio defined on the heavy one.
function _total_averages(curves::Vector{SegmentedCurve}, mass_yields::Vector{MassYield})
    averages = Dict{String,Dict{String,Tuple{Float64,Float64}}}()
    isempty(mass_yields) && return averages
    for curve in curves
        per_distribution = Dict{String,Tuple{Float64,Float64}}()
        for distribution in mass_yields
            try
                per_distribution[distribution.label] = total_average(curve.R_T, distribution)
            catch exception
                exception isa ArgumentError || rethrow()
                @warn "no total average" dataset = curve.label yield = distribution.label reason =
                    exception.msg
            end
        end
        isempty(per_distribution) || (averages[curve.label] = per_distribution)
    end
    return averages
end

# Fit one ratio curve and carry it through to the temperature ratio. Returns nothing when the
# curve cannot support a fit at all, which happens for datasets covering only a few mass pairs.
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
            pinned_value = pinned ? 0.5 : nothing,
            required_windows = windows,
            parsimony = settings.parsimony,
            bounds = (0.0, 1.0),
        )
    catch exception
        exception isa ArgumentError || rethrow()
        @warn "no segmented curve for this dataset" dataset = label reason = exception.msg
        return nothing
    end

    r_ν = evaluate(fit, first(curve.A_H):last(curve.A_H); label = label)
    R_T = temperature_ratio(r_ν, R_a)
    isempty(R_T) && return nothing
    return SegmentedCurve(String(label), fit, r_ν, R_T)
end

# The identities that hold by construction at the symmetric split, checked rather than assumed.
function _symmetry_diagnostics(
    configuration::Configuration,
    domain::FragmentationDomain,
    R_a::AbstractDict{Int,Float64},
    curves::Vector{SegmentedCurve},
)
    warnings = Pair{String,String}[]
    checks = Dict{String,Any}()

    A₀ = configuration.system.A₀
    invariant = symmetric_charge_set_is_invariant(domain, A₀, configuration.system.Z₀)
    checks["symmetric_charge_set_invariant"] =
        invariant === missing ? "not applicable" : invariant
    if invariant === false
        push!(
            warnings,
            "charge_set" => "the charge numbers retained at the symmetric split are not invariant under \
                 Z -> Z₀ - Z, so R_a = 1 and R_T = 1 hold there only approximately; this occurs \
                 when Z₀ is odd and the most probable charge falls between two integers",
        )
    end

    if iseven(A₀)
        A_sym = A₀ ÷ 2
        if haskey(R_a, A_sym)
            deviation = abs(R_a[A_sym] - 1)
            checks["R_a_at_symmetry"] = R_a[A_sym]
            deviation < 1e-8 || push!(
                warnings,
                "R_a" => "R_a departs from unity by $(round(deviation; sigdigits = 3)) at the \
                     symmetric split, where the two fragments are the same nuclide and the ratio \
                     must be exactly one",
            )
        end
        checks["r_nu_pinned_at_symmetry"] = [
            c.label for c in curves if c.fit.pinned_value !== nothing && c.fit.x₀ == A_sym
        ]
    end
    # A pin anywhere but the symmetric split asserts r_ν = 1/2 where no identity holds.
    for c in curves
        c.fit.pinned_value === nothing && continue
        2 * c.fit.x₀ == A₀ || push!(
            warnings,
            "pin" => "the curve \"$(c.label)\" is pinned at A_H = $(c.fit.x₀), which is not the \
                 symmetric split of A₀ = $(A₀)",
        )
    end

    checks["warnings"] = warnings
    return checks
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
    for curve in result.segmented_curves
        label = curve.label
        trend = label == TREND_LABEL
        averages = get(result.total_average_R_T, label, Dict{String,Tuple{Float64,Float64}}())
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
            "range_mean_R_T" => collect(result.range_mean_R_T[label]),
            "total_average_R_T" => Dict{String,Any}(k => collect(v) for (k, v) in averages),
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
        ),
        "segmented_curve" => entries,
    )
    path = _unused_path(joinpath(directory, "manifest_$(identifier).toml"))
    mkpath(dirname(path))
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
# column, so a reader holding the file knows what is in it without consulting the file name.
function _ratio_table(curve::RatioCurve, quantity::AbstractString, significant_digits::Integer)
    return DataFrame(
        :A_H => curve.A_H,
        Symbol(quantity) => round.(curve.ratio; sigdigits = significant_digits),
        Symbol("$(quantity)_uncertainty") => round.(curve.σ; sigdigits = significant_digits),
    )
end

"""
    write_results(result) -> Dict{String,String}

Write the tabulated ratios, the segment pivots, the figures and the run metadata.

Output goes to `<root>/results/<subdirectory>` and `<root>/plots/<subdirectory>`, both named by
the configuration. `root` defaults to the active project, which is what a run from this repository
wants; a caller using the package as a library passes its own. Existing files are never
overwritten; a suffix is appended instead, so that a rerun cannot destroy a previous result.
Returns the paths written, keyed by content.
"""
function write_results(result::ExtractionResult; root::AbstractString = projectdir())
    configuration = result.configuration
    identifier = run_identifier(configuration)
    results_root = joinpath(root, "results", configuration.output.subdirectory)
    plots_root = joinpath(root, "plots", configuration.output.subdirectory)
    digits = configuration.output.significant_digits
    written = Dict{String,String}()

    for (curves, name, quantity) in (
        (result.r_ν, "r_nu_vs_A_H", "r_nu"),
        (result.R_T, "R_T_vs_A_H", "R_T"),
        ([result.consensus_r_ν], "r_nu_vs_A_H_consensus", "r_nu"),
    )
        for curve in curves
            isempty(curve) && continue
            path = _unused_path(
                joinpath(results_root, "$(name)_$(_file_token(curve.label)).csv")
            )
            mkpath(dirname(path))
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
        path = _unused_path(
            joinpath(results_root, "r_nu_vs_A_H_pivots_$(token)_$(identifier).csv")
        )
        mkpath(dirname(path))
        CSV.write(path, pivot_table)
        written["r_nu_pivots/$(curve.label)"] = path

        for (segmented, name, quantity) in (
            (curve.r_ν, "r_nu_vs_A_H_segmented", "r_nu"),
            (curve.R_T, "R_T_vs_A_H_segmented", "R_T"),
        )
            path = _unused_path(joinpath(results_root, "$(name)_$(token)_$(identifier).csv"))
            CSV.write(path, _ratio_table(segmented, quantity, digits))
            written["$(quantity)_segmented/$(curve.label)"] = path
        end
    end

    # The total average over each yield distribution: one row per segmented curve and
    # distribution, which is how the literature tabulates it.
    if !isempty(result.total_average_R_T)
        rows = NamedTuple{
            (:segmented_curve, :mass_yield, :R_T, :R_T_uncertainty),
            Tuple{String,String,Float64,Float64},
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
                        R_T = round(entry[1]; sigdigits = digits),
                        R_T_uncertainty = round(entry[2]; sigdigits = digits),
                    ),
                )
            end
        end
        if !isempty(rows)
            path = _unused_path(joinpath(results_root, "total_average_R_T_$(identifier).csv"))
            mkpath(dirname(path))
            CSV.write(path, DataFrame(rows))
            written["total_average_R_T"] = path
        end
    end

    # Figures come from the CairoMakie extension. A run without it still writes every table and
    # its metadata, rather than failing at the last step for want of a plotting stack.
    if _plotting_extension() === nothing
        @warn "figures skipped; load CairoMakie alongside this package to write them"
    else
        merge!(written, write_figures(result, plots_root, identifier))
    end

    # Per-dataset diagnostics, written for every dataset whether or not it was pooled.
    if !isempty(result.dataset_diagnostics)
        excluded = configuration.excluded_datasets
        rows = [
            (
                dataset = d.label,
                points = d.points,
                pairs = d.pairs,
                first_pair = d.first_pair,
                last_pair = d.last_pair,
                outside_physical_range = d.outside_physical_range,
                symmetry_departure = d.symmetry_departure,
                complement_sum = d.complement_sum,
                complement_spread = d.complement_spread,
                without_uncertainties = d.without_uncertainties,
                pooled = !haskey(excluded, d.label),
                exclusion_reason = get(excluded, d.label, ""),
            ) for d in result.dataset_diagnostics
        ]
        path = _unused_path(joinpath(results_root, "dataset_diagnostics_$(identifier).csv"))
        mkpath(dirname(path))
        CSV.write(path, DataFrame(rows))
        written["dataset_diagnostics"] = path
    end

    written["manifest"] = _write_manifest(result, results_root, identifier, written)

    metadata = run_metadata(configuration)
    metadata["result"] = Dict{String,Any}(
        "segmented_curves" => Dict{String,Any}(
            curve.label => Dict{String,Any}(
                "segments" => segments(curve.fit),
                "breakpoints" => curve.fit.breakpoints,
                "pivots" => [[point[1], point[2]] for point in pivots(curve.fit)],
                "reduced_chi_squared" => curve.fit.wrss / curve.fit.dof,
                "bic" => curve.fit.bic,
                "bic_by_order" =>
                    [[order, value] for (order, value) in curve.fit.selection],
                "weights_imputed" => curve.fit.weights_imputed,
                "range_mean_R_T" => result.range_mean_R_T[curve.label][1],
                "range_mean_R_T_uncertainty" => result.range_mean_R_T[curve.label][2],
            ) for curve in result.segmented_curves
        ),
        "diagnostics" => Dict(k => v for (k, v) in result.diagnostics if k != "warnings"),
    )
    written["metadata"] = write_metadata(
        joinpath(results_root, "metadata_$(identifier).toml"), metadata
    )

    @info "results written" results = results_root plots = plots_root
    return written
end

_file_token(label::AbstractString) = replace(strip(label), r"[^A-Za-z0-9.\-]+" => "_")
