# The pipeline: from experimental multiplicity data to a parameterized temperature ratio.

"""
    PipelineResult

Everything a run produces, held together so that it can be inspected interactively as well as
written to disk.

# Fields

- `configuration`: the configuration the run was driven by.
- `data_sets`: the experimental multiplicity data that was read.
- `r_ν`, `R_T`: the multiplicity and temperature ratios extracted from each data set.
- `R_a`: the level density parameter ratio against heavy-fragment mass number.
- `fit`: the piecewise-linear parameterization of the pooled multiplicity ratio.
- `r_ν_fitted`, `R_T_fitted`: the parameterized multiplicity ratio and the temperature ratio
  obtained from it.
- `range_mean_R_T`: the inverse-variance weighted mean of the parameterized temperature ratio
  over the fragment mass range, with its uncertainty. This is *not* the total average quoted in
  the literature, which is taken over a fission fragment mass yield distribution `Y(A)`;
  computing that requires `Y(A)` as an additional input, which this package does not take.
- `diagnostics`: checks of the exact identities at the symmetric split.
"""
struct PipelineResult
    configuration::Configuration
    data_sets::Vector{MultiplicityData}
    r_ν::Vector{RatioCurve}
    R_T::Vector{RatioCurve}
    R_a::Dict{Int,Float64}
    fit::SegmentedFit
    r_ν_fitted::RatioCurve
    R_T_fitted::RatioCurve
    range_mean_R_T::Tuple{Float64,Float64}
    diagnostics::Dict{String,Any}
end

"""
    read_multiplicity_directory(directory) -> Vector{MultiplicityData}

Read every data file in `directory`, sorted by name so that colours and markers are assigned
reproducibly.

The label of each set is its file name stripped of the leading case identifier and the extension,
with underscores replaced by spaces, which is how the sets are named in the literature.
"""
function read_multiplicity_directory(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("multiplicity directory not found: $(directory)"))
    files = sort!(filter!(f -> !startswith(f, "."), readdir(directory)))
    isempty(files) && throw(ArgumentError("multiplicity directory is empty: $(directory)"))

    data_sets = MultiplicityData[]
    for file in files
        stem = splitext(file)[1]
        label = replace(replace(stem, r"^[A-Za-z0-9]+_[a-z0-9]+f_nuA_" => ""), '_' => ' ')
        # Initials are written without a space in the file names; restore it for the legend.
        label = replace(label, r"(?<=\b[A-Z])\.(?=[A-Z][a-z])" => ". ")
        push!(data_sets, read_multiplicity(joinpath(directory, file); label = label))
    end
    return data_sets
end

"""
    pool(curves; label) -> RatioCurve

Concatenate ratio curves from several data sets into one, sorted by mass number.

The parameterization is fitted to the pooled ratio rather than to each set separately: the sets
of one fissioning nucleus measure the same quantity, and where they disagree that disagreement is
information about the uncertainty of the parameterization, not a reason to choose between them.
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

The steps are: read the mass excesses and the experimental multiplicity data; build the
fragmentation range and the isobaric charge distribution; form the multiplicity ratio
`r_ν = ν_H/(ν_L + ν_H)` of every data set and the temperature ratio from it; fit the pooled
multiplicity ratio by joined straight segments; and obtain the parameterized temperature ratio
from the fitted multiplicity ratio.

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

    pooled = pool(r_ν[usable]; label = "pooled data")
    fit = fit_segments(
        pooled.A_H,
        pooled.value,
        pooled.σ;
        max_segments = configuration.segments.max_segments,
        min_points_per_segment = configuration.segments.min_points_per_segment,
        pinned_value = configuration.segments.pin_symmetric_split ? 0.5 : nothing,
        required_windows = configuration.segments.required_windows,
        bounds = (0.0, 1.0),
    )
    @info "parameterization selected" segments = segments(fit) breakpoints = fit.breakpoints reduced_chi_squared =
        fit.wrss / fit.dof bic = fit.bic

    fitted_range = first(pooled.A_H):last(pooled.A_H)
    r_ν_fitted = evaluate(fit, fitted_range; label = "Segments")
    R_T_fitted = temperature_ratio(r_ν_fitted, R_a)
    isempty(R_T_fitted) && throw(
        ArgumentError("the parameterized multiplicity ratio yields no temperature ratio; the \
                       level density parameter ratio is undefined over the fitted range"),
    )

    diagnostics = _symmetry_diagnostics(configuration, domain, R_a, fit)
    for (key, message) in diagnostics["warnings"]
        @warn message identity = key
    end

    result = PipelineResult(
        configuration,
        data_sets,
        r_ν,
        R_T,
        R_a,
        fit,
        r_ν_fitted,
        R_T_fitted,
        weighted_mean(R_T_fitted),
        diagnostics,
    )

    write_output && write_results(result)
    return result
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

    parameterization = DataFrame(;
        A_H = [pivot[1] for pivot in pivots(result.fit)],
        r_nu = [round(pivot[2]; digits = digits) for pivot in pivots(result.fit)],
    )
    path = _unused_path(joinpath(results_root, "r_nu_segments_$(identifier).csv"))
    mkpath(dirname(path))
    CSV.write(path, parameterization)
    written["r_nu/segments"] = path

    for (curve, name) in
        ((result.r_ν_fitted, "r_nu_parameterized"), (result.R_T_fitted, "R_T_parameterized"))
        table = DataFrame(;
            A_H = curve.A_H,
            value = round.(curve.value; digits = digits),
            uncertainty = round.(curve.σ; digits = digits),
        )
        path = _unused_path(joinpath(results_root, "$(name)_$(identifier).csv"))
        CSV.write(path, table)
        written[name] = path
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
                result.r_ν_fitted;
                ylabel = L"r_\nu = \nu_H / (\nu_L + \nu_H)",
                reference = 0.5,
                reference_label = "Equal sharing",
            ),
        )
        mean_value, mean_uncertainty = result.range_mean_R_T
        return written["figure/R_T"] = save_figure(
            joinpath(plots_root, "R_T_$(identifier).pdf"),
            plot_ratio(
                filter(!isempty, result.R_T),
                result.R_T_fitted;
                ylabel = L"R_T = T_L / T_H",
                reference = 1.0,
                reference_label = "Equal temperatures",
                annotation = "Range mean $(round(mean_value; digits = 3)) ± \
                              $(round(mean_uncertainty; sigdigits = 2))",
            ),
        )
    end

    metadata = run_metadata(configuration)
    metadata["result"] = Dict{String,Any}(
        "segments" => segments(result.fit),
        "breakpoints" => result.fit.breakpoints,
        "reduced_chi_squared" => result.fit.wrss / result.fit.dof,
        "bic" => result.fit.bic,
        "bic_by_order" => [[order, value] for (order, value) in result.fit.selection],
        "weights_imputed" => result.fit.weights_imputed,
        "range_mean_R_T" => result.range_mean_R_T[1],
        "range_mean_R_T_uncertainty" => result.range_mean_R_T[2],
        "diagnostics" => Dict(k => v for (k, v) in result.diagnostics if k != "warnings"),
    )
    written["metadata"] = write_metadata(
        joinpath(results_root, "metadata_$(identifier).toml"), metadata
    )

    @info "results written" results = results_root plots = plots_root
    return written
end

_file_token(label::AbstractString) = replace(strip(label), r"[^A-Za-z0-9.\-]+" => "_")
