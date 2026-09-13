# The pipeline: from experimental multiplicity data to a parameterized temperature ratio.

"""
    Parameterization

One piecewise-linear description of the multiplicity ratio, and the temperature ratio obtained
from it.

A run produces several: one per experimental data set, and one following the systematic behaviour
of the ratio. They are alternatives, not an ensemble to be averaged — a prompt emission code takes
one of them as input, and which one describes reality is settled downstream, by comparing the
multiplicity distributions and yields that code produces against experiment.

# Fields

- `label`: the data set the parameterization came from, or `"systematic trend"`.
- `fit`: the piecewise-linear fit of `r_ν`.
- `r_ν`, `R_T`: the parameterized multiplicity ratio and the temperature ratio from it.
"""
struct Parameterization
    label::String
    fit::SegmentedFit
    r_ν::RatioCurve
    R_T::RatioCurve
end

"""
    PipelineResult

Everything a run produces, held together so that it can be inspected interactively as well as
written to disk.

# Fields

- `configuration`: the configuration the run was driven by.
- `data_sets`: the experimental multiplicity data that was read.
- `r_ν`, `R_T`: the multiplicity and temperature ratios extracted from each data set.
- `R_a`: the level density parameter ratio against heavy-fragment mass number.
- `parameterizations`: one per data set that supports a fit, followed by the systematic-trend
  curve. Alternatives offered to a prompt emission code, not an ensemble.
- `range_mean_R_T`: for each parameterization, the inverse-variance weighted mean of its
  temperature ratio over the fragment mass range, with its uncertainty. This is *not* the total
  average quoted in the literature: it weights every mass number equally and so is dominated by
  the far-asymmetric tail, where the yield is negligible.
- `yield_sets`: the fragment mass yield distributions read, empty when the configuration names
  none.
- `total_average_R_T`: the quantity the literature quotes, `⟨R_T⟩` over each yield distribution,
  keyed by parameterization label and then by yield label. Empty when no yields were given.
- `diagnostics`: checks of the exact identities at the symmetric split.
"""
struct PipelineResult
    configuration::Configuration
    data_sets::Vector{MultiplicityData}
    r_ν::Vector{RatioCurve}
    R_T::Vector{RatioCurve}
    R_a::Dict{Int,Float64}
    parameterizations::Vector{Parameterization}
    range_mean_R_T::Dict{String,Tuple{Float64,Float64}}
    yield_sets::Vector{YieldData}
    total_average_R_T::Dict{String,Dict{String,Tuple{Float64,Float64}}}
    diagnostics::Dict{String,Any}
end

"""
    systematic_trend(result) -> Parameterization

The parameterization that follows the systematic behaviour of the ratio rather than any single
data set.
"""
function systematic_trend(result::PipelineResult)
    index = findfirst(p -> p.label == TREND_LABEL, result.parameterizations)
    index === nothing && throw(ArgumentError("this result carries no systematic-trend curve"))
    return result.parameterizations[index]
end

"""
    read_multiplicity_directory(directory) -> Vector{MultiplicityData}

Read every `.dat` file in `directory`, sorted by name so that colours and markers are assigned
reproducibly. Other files are ignored, so a run record can sit beside the data it describes.

The label of each set is its file name stripped of the leading archive identifier and the
extension, with underscores replaced by spaces, which is how the sets are named in the literature.
"""
function read_multiplicity_directory(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("multiplicity directory not found: $(directory)"))
    files = _data_files(directory)
    isempty(files) &&
        throw(ArgumentError("multiplicity directory holds no data files: $(directory)"))
    return [
        read_multiplicity(joinpath(directory, file); label = _set_label(file)) for file in files
    ]
end

# Data files carry the `.dat` extension; anything else in the directory — a run record, a note —
# is not data and is skipped rather than parsed and rejected.
function _data_files(directory::AbstractString)
    return sort!(filter!(f -> endswith(f, ".dat") && !startswith(f, "."), readdir(directory)))
end

# `<identifier>_<author>_<year>.dat`, or the older `<case>_nuA_<author>_<year>.dat`.
function _set_label(file::AbstractString)
    stem = splitext(file)[1]
    stem = replace(stem, r"^[0-9]+_" => "", r"^[A-Za-z0-9]+_[a-z0-9]+f_nuA_" => "")
    label = replace(stem, '_' => ' ')
    # Initials are written without a space in the file names; restore it for the legend.
    return replace(label, r"(?<=\b[A-Z])\.(?=[A-Z][a-z])" => ". ")
end

"""
    pool(curves; label) -> RatioCurve

Concatenate ratio curves from several data sets into one, sorted by mass number.

Used for the systematic-trend curve, which is meant to follow the general behaviour of the ratio
rather than any one measurement. The per-set parameterizations are fitted to their own data: where
the sets of one fissioning nucleus disagree, as they do markedly for some nuclei, a curve fitted
through all of them at once describes none of them, and a prompt emission code has nothing to
validate it against.
"""
function pool(curves::Vector{RatioCurve}; label::AbstractString = "pooled")
    A_H = reduce(vcat, (curve.A_H for curve in curves); init = Int[])
    value = reduce(vcat, (curve.value for curve in curves); init = Float64[])
    σ = reduce(vcat, (curve.σ for curve in curves); init = Float64[])
    order = sortperm(A_H)
    return RatioCurve(A_H[order], value[order], σ[order], String(label))
end

"""
    run_pipeline(configuration; write_output = true) -> PipelineResult

Execute the extraction for one configuration.

The steps are: read the tabulated data and the experimental multiplicity data; build the
fragmentation range and the isobaric charge distribution; form the multiplicity ratio
`r_ν = ν_H/(ν_L + ν_H)` of every data set and the temperature ratio from it; and parameterize the
multiplicity ratio by joined straight segments, once per data set and once more following its
systematic behaviour.

A run therefore returns several parameterizations rather than one. Where the data sets of a
fissioning nucleus disagree — as they do markedly for some — a curve fitted through all of them
describes none of them, and a prompt emission code has nothing to validate it against. They are
offered as alternatives: the code takes one as input, and which one describes reality is settled
by comparing the multiplicity distributions and yields it produces against experiment. The
systematic-trend curve is the one to use where a data set is too sparse or too scattered to
determine the shape on its own.

The temperature ratio is parameterized by transforming the fitted multiplicity ratio, not by
fitting segments to the temperature ratio a second time. The transformation is exact, whereas a
second fit would discard the uncertainty of the first and impose a piecewise-linear shape on a
quantity that is not piecewise-linear.

With `write_output`, tabulated results, figures and run metadata are written under the results
and plots directories; uncertainty of the input is preserved by never overwriting an existing
file.
"""
function run_pipeline(configuration::Configuration; write_output::Bool = true)
    @info "reading input" configuration = configuration.source
    prescription = build_prescription(configuration.level_density)
    charge_data = read_charge_distribution(
        configuration.fragmentation.charge_distribution_file,
        configuration.fragmentation.fallback_charge_polarization,
        configuration.fragmentation.fallback_rms,
    )
    data_sets = read_multiplicity_directory(configuration.multiplicity_directory)

    A₀ = configuration.system.A₀
    Z₀ = configuration.system.Z₀
    range = A_H_range(configuration)

    domain = fragmentation_domain(
        A₀, Z₀, range, configuration.fragmentation.charges_per_mass, charge_data
    )
    if !isempty(domain.fallback_masses)
        @warn "tabulated charge polarization and dispersion unavailable; the average values were \
               used" masses = _summarize_range(domain.fallback_masses) ΔZ =
            configuration.fragmentation.fallback_charge_polarization rms =
            configuration.fragmentation.fallback_rms
    end

    R_a = level_density_ratio(
        configuration.level_density.ratio_averaging, prescription, A₀, Z₀, domain
    )
    isempty(R_a) &&
        throw(ArgumentError("the level density parameter ratio is undefined over the whole \
                       fragmentation range; check the mass excess table"))

    r_ν = [multiplicity_ratio(data, A₀, range) for data in data_sets]
    R_T = [temperature_ratio(curve, R_a) for curve in r_ν]
    usable = findall(!isempty, r_ν)
    isempty(usable) &&
        throw(ArgumentError("no data set provides both fragments of any pair within \
                       $(first(range)):$(last(range))"))

    settings = configuration.segments
    parameterizations = Parameterization[]

    # One curve per data set. These are the alternatives a prompt emission code chooses between.
    for index in usable
        curve = r_ν[index]
        parameterization = _parameterize(curve, R_a, settings, curve.label, UnitRange{Int}[])
        parameterization === nothing && continue
        push!(parameterizations, parameterization)
    end

    # One curve following the systematic behaviour of the ratio: the minimum at the heavy magic
    # fragment placed rather than fitted, and the rise above the most probable fragmentation taken
    # through the whole body of data, so its slope falls between those of the individual sets.
    # This is the curve to use where a data set is too sparse or too scattered to determine the
    # shape on its own.
    pooled = pool(r_ν[usable]; label = TREND_LABEL)
    trend = _parameterize(pooled, R_a, settings, TREND_LABEL, settings.required_windows)
    if trend === nothing && !isempty(settings.required_windows)
        @warn "no parameterization satisfies the required windows; the systematic-trend curve \
               was fitted without them" windows = settings.required_windows
        trend = _parameterize(pooled, R_a, settings, TREND_LABEL, UnitRange{Int}[])
    end
    trend === nothing || push!(parameterizations, trend)

    isempty(parameterizations) &&
        throw(ArgumentError("no data set supports a parameterization with \
                       min_points_per_segment = $(settings.min_points_per_segment)"))

    for parameterization in parameterizations
        @info "parameterization" set = parameterization.label segments = segments(
            parameterization.fit
        ) breakpoints = parameterization.fit.breakpoints reduced_chi_squared =
            parameterization.fit.wrss / parameterization.fit.dof
    end

    diagnostics = _symmetry_diagnostics(
        configuration, domain, R_a, first(parameterizations).fit
    )
    for (key, message) in diagnostics["warnings"]
        @warn message identity = key
    end

    yield_sets = if configuration.yield_directory === nothing
        YieldData[]
    else
        read_yield_directory(configuration.yield_directory)
    end
    total_average_R_T = _total_averages(parameterizations, yield_sets)
    for parameterization in parameterizations
        for distribution in yield_sets
            entry = get(
                get(total_average_R_T, parameterization.label, Dict()),
                distribution.label,
                nothing,
            )
            entry === nothing && continue
            @info "total average" set = parameterization.label yield = distribution.label R_T = entry[1] uncertainty = entry[2]
        end
    end

    result = PipelineResult(
        configuration,
        data_sets,
        r_ν,
        R_T,
        R_a,
        parameterizations,
        Dict(p.label => weighted_mean(p.R_T) for p in parameterizations),
        yield_sets,
        total_average_R_T,
        diagnostics,
    )

    write_output && write_results(result)
    return result
end

# The total average of every parameterization over every yield distribution. A pair that shares no
# mass number is omitted rather than reported as zero: a distribution covering only the light wing
# says nothing about a ratio defined on the heavy one.
function _total_averages(
    parameterizations::Vector{Parameterization}, yield_sets::Vector{YieldData}
)
    averages = Dict{String,Dict{String,Tuple{Float64,Float64}}}()
    isempty(yield_sets) && return averages
    for parameterization in parameterizations
        per_distribution = Dict{String,Tuple{Float64,Float64}}()
        for distribution in yield_sets
            try
                per_distribution[distribution.label] = total_average(
                    parameterization.R_T, distribution
                )
            catch err
                err isa ArgumentError || rethrow()
                @warn "no total average" set = parameterization.label yield = distribution.label reason =
                    err.msg
            end
        end
        isempty(per_distribution) || (averages[parameterization.label] = per_distribution)
    end
    return averages
end

# Fit one ratio curve and carry it through to the temperature ratio. Returns nothing when the
# curve cannot support a fit at all, which happens for data sets covering only a few mass pairs.
function _parameterize(
    curve::RatioCurve,
    R_a::AbstractDict{Int,Float64},
    settings::SegmentSettings,
    label::AbstractString,
    windows::Vector{UnitRange{Int}},
)
    fit = try
        fit_segments(
            curve.A_H,
            curve.value,
            curve.σ;
            max_segments = settings.max_segments,
            min_points_per_segment = settings.min_points_per_segment,
            pinned_value = settings.pin_symmetric_split ? 0.5 : nothing,
            required_windows = windows,
            bounds = (0.0, 1.0),
        )
    catch err
        err isa ArgumentError || rethrow()
        @warn "no parameterization for this data set" set = label reason = err.msg
        return nothing
    end

    r_ν = evaluate(fit, first(curve.A_H):last(curve.A_H); label = label)
    R_T = temperature_ratio(r_ν, R_a)
    isempty(R_T) && return nothing
    return Parameterization(String(label), fit, r_ν, R_T)
end

# The identities that hold by construction at the symmetric split, checked rather than assumed.
function _symmetry_diagnostics(
    configuration::Configuration,
    domain::FragmentationDomain,
    R_a::AbstractDict{Int,Float64},
    fit::SegmentedFit,
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
        if first((x -> [fit.x₀; x])(fit.breakpoints)) == A_sym && fit.pinned_value !== nothing
            checks["r_nu_pinned_at_symmetry"] = fit.pinned_value
        end
    end

    checks["warnings"] = warnings
    return checks
end

function _summarize_range(masses::Vector{Int})
    isempty(masses) && return "none"
    return "$(first(masses)):$(last(masses)) ($(length(masses)) mass numbers)"
end

"""
    write_results(result) -> Dict{String,String}

Write the tabulated ratios, the segment parameterization, the figures and the run metadata.

Output goes to `results/<subdirectory>` and `plots/<subdirectory>`, both named by the
configuration. Existing files are never overwritten; a suffix is appended instead, so that a
rerun cannot destroy a previous result. Returns the paths written, keyed by content.
"""
function write_results(result::PipelineResult)
    configuration = result.configuration
    identifier = run_identifier(configuration)
    results_root = joinpath(projectdir(), "results", configuration.output.subdirectory)
    plots_root = joinpath(projectdir(), "plots", configuration.output.subdirectory)
    digits = configuration.output.digits
    written = Dict{String,String}()

    for (curves, name) in ((result.r_ν, "r_nu"), (result.R_T, "R_T"))
        for curve in curves
            isempty(curve) && continue
            table = DataFrame(;
                A_H = curve.A_H,
                value = round.(curve.value; digits = digits),
                uncertainty = round.(curve.σ; digits = digits),
            )
            path = _unused_path(
                joinpath(results_root, "$(name)_$(_file_token(curve.label)).csv")
            )
            mkpath(dirname(path))
            CSV.write(path, table)
            written["$(name)/$(curve.label)"] = path
        end
    end

    for parameterization in result.parameterizations
        token = _file_token(parameterization.label)
        points = pivots(parameterization.fit)
        segment_table = DataFrame(;
            A_H = [point[1] for point in points],
            r_nu = [round(point[2]; digits = digits) for point in points],
        )
        path = _unused_path(joinpath(results_root, "segments_$(token)_$(identifier).csv"))
        mkpath(dirname(path))
        CSV.write(path, segment_table)
        written["segments/$(parameterization.label)"] = path

        for (curve, name) in (
            (parameterization.r_ν, "r_nu_parameterized"),
            (parameterization.R_T, "R_T_parameterized"),
        )
            table = DataFrame(;
                A_H = curve.A_H,
                value = round.(curve.value; digits = digits),
                uncertainty = round.(curve.σ; digits = digits),
            )
            path = _unused_path(joinpath(results_root, "$(name)_$(token)_$(identifier).csv"))
            CSV.write(path, table)
            written["$(name)/$(parameterization.label)"] = path
        end
    end

    # The total average over each yield distribution: one row per parameterization and
    # distribution, which is how the literature tabulates it.
    if !isempty(result.total_average_R_T)
        rows = NamedTuple{
            (:parameterization, :yield_distribution, :R_T, :uncertainty),
            Tuple{String,String,Float64,Float64},
        }[]
        for parameterization in result.parameterizations
            per_distribution = get(result.total_average_R_T, parameterization.label, nothing)
            per_distribution === nothing && continue
            for distribution in result.yield_sets
                entry = get(per_distribution, distribution.label, nothing)
                entry === nothing && continue
                push!(
                    rows,
                    (
                        parameterization = parameterization.label,
                        yield_distribution = distribution.label,
                        R_T = round(entry[1]; digits = digits),
                        uncertainty = round(entry[2]; digits = digits),
                    ),
                )
            end
        end
        if !isempty(rows)
            path = _unused_path(joinpath(results_root, "total_average_$(identifier).csv"))
            mkpath(dirname(path))
            CSV.write(path, DataFrame(rows))
            written["total_average"] = path
        end
    end

    with_theme(publication_theme()) do
        written["figure/multiplicity"] = save_figure(
            joinpath(plots_root, "multiplicity_$(identifier).pdf"),
            plot_multiplicities(result.data_sets; A₀ = configuration.system.A₀),
        )
        written["figure/r_nu"] = save_figure(
            joinpath(plots_root, "r_nu_$(identifier).pdf"),
            plot_ratio(
                filter(!isempty, result.r_ν),
                [p.r_ν for p in result.parameterizations];
                ylabel = L"r_\nu = \nu_H / (\nu_L + \nu_H)",
                reference = 0.5,
                reference_label = "Equal sharing",
            ),
        )
        return written["figure/R_T"] = save_figure(
            joinpath(plots_root, "R_T_$(identifier).pdf"),
            plot_ratio(
                filter(!isempty, result.R_T),
                [p.R_T for p in result.parameterizations];
                ylabel = L"R_T = T_L / T_H",
                reference = 1.0,
                reference_label = "Equal temperatures",
            ),
        )
    end

    metadata = run_metadata(configuration)
    metadata["result"] = Dict{String,Any}(
        "parameterizations" => Dict{String,Any}(
            parameterization.label => Dict{String,Any}(
                "segments" => segments(parameterization.fit),
                "breakpoints" => parameterization.fit.breakpoints,
                "pivots" =>
                    [[point[1], point[2]] for point in pivots(parameterization.fit)],
                "reduced_chi_squared" =>
                    parameterization.fit.wrss / parameterization.fit.dof,
                "bic" => parameterization.fit.bic,
                "bic_by_order" => [
                    [order, value] for (order, value) in parameterization.fit.selection
                ],
                "weights_imputed" => parameterization.fit.weights_imputed,
                "range_mean_R_T" => result.range_mean_R_T[parameterization.label][1],
                "range_mean_R_T_uncertainty" =>
                    result.range_mean_R_T[parameterization.label][2],
            ) for parameterization in result.parameterizations
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
