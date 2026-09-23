# Figures in the standard layout: a 900 × 600 canvas per panel, widened to 1200 when many datasets
# share an axis, 26 pt type, 3 pt data lines and 14 pt stroked markers.
#
# The implementations live here rather than in the package so that `using FissionTemperatureRatio`
# does not load a plotting stack. CairoMakie depends on LaTeXStrings and MathTeXEngine, so loading
# CairoMakie alone triggers this extension. The docstrings are on the stubs in `src/plotting.jl`,
# where the documentation build can find them.

module FissionTemperatureRatioCairoMakieExt

using CairoMakie
using Statistics: quantile
using LaTeXStrings: @L_str, LaTeXString
using Printf: @sprintf
using MathTeXEngine: texfont

using FissionTemperatureRatio:
    FissionTemperatureRatio, ExtractionResult, Multiplicity, RatioCurve, TREND_LABEL

const BASE_WIDTH = 900
const WIDE_WIDTH = 1200
const PANEL_HEIGHT = 600
const LEGEND_ROW_HEIGHT = 38
# Above this many legend entries the canvas widens, so the legend banks into three columns.
const WIDE_ABOVE_ENTRIES = 8
const ANNOTATION_SIZE = 21

function FissionTemperatureRatio.publication_theme()
    return Theme(;
        fonts = (; regular = texfont(:text), bold = texfont(:bold), italic = texfont(:italic)),
        fontsize = 26,
        figure_padding = 10,
        linewidth = 3,
        markersize = 14,
        Axis = (
            spinewidth = 1.5,
            xtickwidth = 1.5,
            ytickwidth = 1.5,
            xticklabelsize = 22,
            yticklabelsize = 22,
            xlabelpadding = 8,
            ylabelpadding = 8,
            xgridstyle = :dash,
            ygridstyle = :dash,
            xgridcolor = (:grey, 0.12),
            ygridcolor = (:grey, 0.12),
            xminorticksvisible = false,
            yminorticksvisible = false,
            xtickalign = 1,
            ytickalign = 1,
        ),
        Scatter = (strokewidth = 1.5,),
        Legend = (framevisible = false, orientation = :horizontal, titlefont = :bold),
    )
end

# One colour per dataset, consistent across every figure of a run. Okabe-Ito, which stays
# distinguishable in grayscale and for the common colour vision deficiencies.
const DATASET_COLORS = [
    RGBf(0.0, 0.447, 0.698),
    RGBf(0.835, 0.369, 0.0),
    RGBf(0.0, 0.620, 0.451),
    RGBf(0.800, 0.475, 0.655),
    RGBf(0.337, 0.706, 0.914),
    RGBf(0.941, 0.894, 0.259),
    RGBf(0.902, 0.624, 0.0),
    RGBf(0.35, 0.35, 0.35),
]

# Seven markers against eight colours, deliberately: the two cycles are then coprime, so the
# (colour, marker) pair is unique for 56 datasets rather than repeating every eight. With sixteen
# datasets for one nucleus, equal cycles put two measurements under the same blue circle.
const DATASET_MARKERS = [:circle, :rect, :utriangle, :diamond, :dtriangle, :xcross, :star5]

dataset_color(index::Integer) = DATASET_COLORS[mod1(index, length(DATASET_COLORS))]
dataset_marker(index::Integer) = DATASET_MARKERS[mod1(index, length(DATASET_MARKERS))]

# A dataset keeps one colour and one marker across every figure of a run. Indexing by position
# cannot do that: the multiplicity figure is handed every dataset, the ratio figures only those
# with a usable ratio, and the fitted curves only those that supported a fit — three lists of
# different lengths, so position means a different measurement in each. One dataset failing to fit
# shifted every later curve onto its neighbour's colour. The order of the labels as read is the
# key instead, and an unknown label falls back to its position.
function _style_index(label::AbstractString, order::Vector{String}, fallback::Integer)
    index = findfirst(==(String(label)), order)
    return index === nothing ? fallback : index
end

figure_width(entries::Integer) = entries > WIDE_ABOVE_ENTRIES ? WIDE_WIDTH : BASE_WIDTH
# Author-and-year labels at 26 pt: two fit across the base canvas, three across the wide one.
legend_columns(entries::Integer) = entries > WIDE_ABOVE_ENTRIES ? 3 : 2
legend_rows(entries::Integer) = max(1, cld(entries, legend_columns(entries)))

# The axes keep a constant height; the figure grows to make room for the legend above them.
function _figure_size(entries::Integer)
    return (figure_width(entries), PANEL_HEIGHT + LEGEND_ROW_HEIGHT * legend_rows(entries))
end

# Marker edge: the fill colour darkened, so markers stay separable where they overlap.
_darker(color::RGBf) = RGBf(0.6f0 * color.r, 0.6f0 * color.g, 0.6f0 * color.b)

# Error bars on the points that quote an uncertainty; a point quoting none is drawn bare.
function _errorbars!(axis::Axis, x, y, σ, color)
    quoted = findall(!ismissing, σ)
    isempty(quoted) && return nothing
    errorbars!(
        axis,
        x[quoted],
        y[quoted],
        Float64[σ[i] for i in quoted];
        color = color,
        linewidth = 1.5,
        whiskerwidth = 0,
    )
    return nothing
end

# The band of a fitted curve, whose uncertainties are all quoted: zero at a pin, positive
# elsewhere.
_band_width(curve::RatioCurve) = Float64[coalesce(σ, 0.0) for σ in curve.σ]

# Mass-number ticks at multiples of ten inside the limits.
_mass_ticks(low::Real, high::Real) = (10 * cld(floor(Int, low), 10)):10:floor(Int, high)

function FissionTemperatureRatio.plot_multiplicities(
    datasets::Vector{Multiplicity}; A_0::Integer, order::Vector{String} = String[]
)
    figure = Figure(; size = _figure_size(length(datasets)))
    axis = Axis(
        figure[1, 1];
        xlabel = L"Fragment mass number $A$",
        ylabel = L"Prompt neutron multiplicity $\nu$",
    )

    for (position, data) in enumerate(datasets)
        index = _style_index(data.label, order, position)
        color = dataset_color(index)
        _errorbars!(axis, data.A, data.ν, data.σν, color)
        scatter!(
            axis,
            data.A,
            data.ν;
            color = color,
            strokecolor = _darker(color),
            marker = dataset_marker(index),
            label = data.label,
        )
    end

    # The far-asymmetric tail of some measurements reaches tens of neutrons per fragment with
    # uncertainties to match. Scaled to those, the sawtooth that carries the physics collapses to a
    # flat line, so the view is bounded by the bulk of the data. No point is discarded — points
    # above the bound simply fall outside the axes.
    bulk = reduce(vcat, (data.ν for data in datasets); init = Float64[])
    if !isempty(bulk)
        ylims!(axis, 0, 1.15 * quantile(bulk, 0.99))
    end
    isempty(datasets) || (
        axis.xticks = _mass_ticks(
            minimum(minimum(data.A) for data in datasets),
            maximum(maximum(data.A) for data in datasets),
        )
    )

    _legend_above(figure, axis, length(datasets))
    return figure
end

function FissionTemperatureRatio.plot_ratio(
    curves::Vector{RatioCurve},
    fitted::Vector{RatioCurve} = RatioCurve[];
    ylabel,
    reference::Union{Real,Nothing} = nothing,
    reference_label::AbstractString = "",
    order::Vector{String} = String[],
    limits = nothing,
    annotation::Union{AbstractString,AbstractVector{<:AbstractString}} = "",
    fit_label::AbstractString = "",
    yticks = nothing,
    annotation_corner::Symbol = :rt,
)
    # A fitted curve reaches the legend if it is the trend, which is always named, or if the
    # caller named the dataset's own fit.
    fitted_labelled =
        curve -> !isempty(curve) && (curve.label == TREND_LABEL || !isempty(fit_label))
    entries =
        count(!isempty, curves) +
        (reference !== nothing && !isempty(reference_label) ? 1 : 0) +
        count(fitted_labelled, fitted)
    # Alone, the trend carries the figure; beside a dataset's own fit it is a guide behind it.
    trend_width = count(!isempty, fitted) == 1 ? 3 : 1.5

    figure = Figure(; size = _figure_size(entries))
    axis = Axis(figure[1, 1]; xlabel = L"Heavy fragment mass number $A_H$", ylabel = ylabel)
    # Linked panels of one run share an abscissa, and a ratio bounded by construction is shown
    # against its bounds rather than against whatever the data happened to reach.
    limits === nothing || limits!(axis, limits...)

    if reference !== nothing
        hlines!(
            axis,
            [reference];
            color = :black,
            linestyle = :dash,
            linewidth = 1.5,
            label = isempty(reference_label) ? nothing : reference_label,
        )
    end

    # The bands go down first, all of them, so that no band hides another curve's markers.
    for (position, curve) in enumerate(fitted)
        isempty(curve) && continue
        trend = curve.label == TREND_LABEL
        color =
            trend ? RGBf(0, 0, 0) : dataset_color(_style_index(curve.label, order, position))
        # Both ratios drawn here are non-negative by construction, so the band is clipped at zero
        # rather than drawn into a region the quantity cannot occupy. The symmetric interval is a
        # Gaussian approximation; where it reaches below zero it is the approximation failing, not
        # the quantity.
        width = _band_width(curve)
        band!(
            axis,
            curve.A_H,
            max.(curve.ratio .- width, 0.0),
            curve.ratio .+ width;
            color = (color, trend ? 0.25 : 0.35),
        )
    end

    for (position, curve) in enumerate(curves)
        isempty(curve) && continue
        index = _style_index(curve.label, order, position)
        color = dataset_color(index)
        _errorbars!(axis, curve.A_H, curve.ratio, curve.σ, color)
        scatter!(
            axis,
            curve.A_H,
            curve.ratio;
            color = color,
            strokecolor = _darker(color),
            marker = dataset_marker(index),
            label = curve.label,
        )
    end

    # The segmented curves are alternatives, so each is drawn in the colour of the dataset it came
    # from, and the systematic-trend curve in black, dashed, to mark that it follows no single
    # measurement.
    for (position, curve) in enumerate(fitted)
        isempty(curve) && continue
        trend = curve.label == TREND_LABEL
        color =
            trend ? RGBf(0, 0, 0) : dataset_color(_style_index(curve.label, order, position))
        label = if trend
            uppercasefirst(TREND_LABEL)
        else
            isempty(fit_label) ? nothing : fit_label
        end
        lines!(
            axis,
            curve.A_H,
            curve.ratio;
            color = color,
            linewidth = trend ? trend_width : 3,
            linestyle = trend ? :dash : :solid,
            label = label,
        )
    end

    # The takeaway belongs where the reader is looking, inside the axes.
    _annotate!(axis, annotation; corner = annotation_corner)

    # Mass numbers read best at decades; the bounds are the ones the axis actually shows.
    bounds = if limits !== nothing && limits[1] isa Real && limits[2] isa Real
        (limits[1], limits[2])
    else
        drawn = filter(!isempty, curves)
        if isempty(drawn)
            nothing
        else
            (minimum(minimum(c.A_H) for c in drawn), maximum(maximum(c.A_H) for c in drawn))
        end
    end
    bounds === nothing || (axis.xticks = _mass_ticks(bounds...))
    yticks === nothing || (axis.yticks = yticks)

    _legend_above(figure, axis, entries)
    return figure
end

# An in-axis annotation, one string per line. Each line is drawn on its own because a line is one
# typeset expression and MathTeXEngine has no line break; the leading is therefore set in the same
# units as the font, from a shared anchor, which a relative coordinate cannot express.
#
# The anchor is the upper right, the one corner a temperature ratio leaves free: R_T is unity at
# the symmetric split, peaks near the shell closure and falls below unity towards the heavy wing,
# so the large-mass, large-ratio corner holds no data in any system. The lower right, where this
# annotation used to sit, is exactly where the heavy wing descends through it.
function _annotate!(
    axis::Axis, annotation; corner::Symbol = :rt, fontsize = ANNOTATION_SIZE, leading = 1.3
)
    corner in (:rt, :rb) ||
        throw(ArgumentError("annotation corner must be :rt or :rb, got $(repr(corner))"))
    lines =
        annotation isa AbstractString ? (isempty(annotation) ? () : (annotation,)) : annotation
    top = corner === :rt
    for (position, line) in enumerate(lines)
        # Stacked downwards from the top corner, upwards from the bottom one.
        shift = top ? -(position - 1) : length(lines) - position
        text!(
            axis,
            0.98,
            top ? 0.96 : 0.04;
            text = line,
            space = :relative,
            align = (:right, top ? :top : :bottom),
            offset = (0, shift * leading * fontsize),
            fontsize = fontsize,
        )
    end
    return axis
end

# Legends sit above the axes, horizontally, so that they cannot collide with the data however the
# points and their error bars happen to fall.
function _legend_above(figure::Figure, axis::Axis, entries::Integer)
    # Banked so that no row exceeds the column count the canvas can hold; wider rows overflow the
    # figure and the rightmost label is clipped without any error.
    Legend(
        figure[0, 1],
        axis;
        orientation = :horizontal,
        nbanks = legend_rows(entries),
        framevisible = false,
        padding = (0, 0, 0, 0),
        tellheight = true,
        tellwidth = false,
    )
    rowgap!(figure.layout, 10)
    return figure
end

function FissionTemperatureRatio.save_figure(path::AbstractString, figure::Figure)
    mkpath(dirname(path))
    target = FissionTemperatureRatio._unused_path(path)
    save(target, figure)
    return target
end

# A value and its uncertainty share their decimal places, as a measured quantity is written: three
# of them, which is what a ratio near unity with an uncertainty of a few parts per thousand needs.
# `round` alone would drop a trailing zero and print 1.18 ± 0.006, two precisions for one number.
_measured(value::Real, uncertainty::Real) = @sprintf("%.3f \\pm %.3f", value, uncertainty)

# The systematic-trend result, as the reader of the figure wants it: the total average where a
# yield distribution was given, and the range mean otherwise, since only one of them exists.
function _trend_annotation(result::ExtractionResult)
    haskey(result.range_mean_R_T, TREND_LABEL) || return LaTeXString[]
    averages = get(result.total_average_R_T, TREND_LABEL, nothing)
    if averages !== nothing && !isempty(averages)
        name = first(sort!(collect(keys(averages))))
        value, uncertainty = averages[name].value, averages[name].uncertainty
        return [
            L"Trend $\langle R_T \rangle = %$(_measured(value, uncertainty))$",
            L"over $Y(A)$ of %$(name)",
        ]
    end
    value, uncertainty = result.range_mean_R_T[TREND_LABEL]
    return [L"Trend range mean $= %$(_measured(value, uncertainty))$"]
end

# Order and quality of one segmented fit, for the figure of that dataset. A dataset quoting no
# uncertainties has uniform weights, and its residual sum is not a chi-squared; the figure says so
# rather than printing a number that means nothing.
function _fit_annotation(curve)
    fit = curve.fit
    count = FissionTemperatureRatio.segments(fit)
    points = fit.dof + length(fit.coefficients) + length(fit.breakpoints)
    quality = if fit.weights_imputed == points
        LaTeXString("no quoted uncertainties")
    else
        reduced = @sprintf("%.2f", fit.wrss / fit.dof)
        L"$\chi^2/\mathrm{dof} = %$(reduced)$"
    end
    return [LaTeXString(count == 1 ? "1 segment" : "$(count) segments"), quality]
end

# The run's figures. Held here, rather than in the pipeline, so that the pipeline carries no
# reference to a plotting type; `write_results` calls it through `Base.get_extension`.
#
# Three overview figures, and then one pair per dataset that supports a fit: on a single axis the
# overview cannot show a dozen fitted curves without burying the measurements under them, so each
# fit is shown against its own data, with the systematic trend beside it for comparison.
function FissionTemperatureRatio.write_figures(
    result::ExtractionResult, directory::AbstractString
)
    written = Dict{String,String}()
    configuration = result.configuration
    order = [data.label for data in result.datasets]
    masses = FissionTemperatureRatio.A_H_range(configuration)
    trend_r_ν = [c.r_ν for c in result.segmented_curves if c.label == TREND_LABEL]
    trend_R_T = [c.R_T for c in result.segmented_curves if c.label == TREND_LABEL]

    # The temperature ratio is unbounded above and its far-asymmetric tail runs away, so the view
    # is bounded by the bulk of the measurements, as the multiplicity figure is.
    R_T_data = filter(!isempty, result.R_T)
    R_T_values = reduce(vcat, (curve.ratio for curve in R_T_data); init = Float64[])
    R_T_upper = isempty(R_T_values) ? nothing : 1.15 * quantile(R_T_values, 0.99)
    R_T_lower = R_T_upper === nothing ? nothing : 0

    with_theme(FissionTemperatureRatio.publication_theme()) do
        written["figure/nu_vs_A"] = FissionTemperatureRatio.save_figure(
            joinpath(directory, "nu_vs_A.pdf"),
            FissionTemperatureRatio.plot_multiplicities(
                result.datasets; A_0 = configuration.system.A₀, order = order
            ),
        )
        written["figure/r_nu_vs_A_H"] = FissionTemperatureRatio.save_figure(
            joinpath(directory, "r_nu_vs_A_H.pdf"),
            FissionTemperatureRatio.plot_ratio(
                filter(!isempty, result.r_ν),
                trend_r_ν;
                ylabel = L"r_\nu = \nu_H / (\nu_L + \nu_H)",
                reference = 0.5,
                reference_label = "Equal sharing",
                order = order,
                limits = (first(masses) - 1, last(masses) + 1, 0, 1),
                yticks = 0:0.25:1,
            ),
        )
        written["figure/R_T_vs_A_H"] = FissionTemperatureRatio.save_figure(
            joinpath(directory, "R_T_vs_A_H.pdf"),
            FissionTemperatureRatio.plot_ratio(
                R_T_data,
                trend_R_T;
                ylabel = L"R_T = T_L / T_H",
                reference = 1.0,
                reference_label = "Equal temperatures",
                order = order,
                limits = (first(masses) - 1, last(masses) + 1, R_T_lower, R_T_upper),
                annotation = _trend_annotation(result),
            ),
        )

        for curve in result.segmented_curves
            curve.label == TREND_LABEL && continue
            token = FissionTemperatureRatio._file_token(curve.label)
            annotation = _fit_annotation(curve)

            index = findfirst(c -> c.label == curve.label, result.r_ν)
            if index !== nothing && !isempty(result.r_ν[index])
                written["figure/r_nu_vs_A_H/$(curve.label)"] = FissionTemperatureRatio.save_figure(
                    joinpath(directory, "r_nu_vs_A_H_segmented_$(token).pdf"),
                    FissionTemperatureRatio.plot_ratio(
                        RatioCurve[result.r_ν[index]],
                        vcat(RatioCurve[curve.r_ν], trend_r_ν);
                        ylabel = L"r_\nu = \nu_H / (\nu_L + \nu_H)",
                        reference = 0.5,
                        reference_label = "Equal sharing",
                        order = order,
                        limits = (first(masses) - 1, last(masses) + 1, 0, 1),
                        annotation = annotation,
                        # r_ν rises into the upper right; the lower right is the free corner.
                        annotation_corner = :rb,
                        fit_label = "Segmented fit",
                        yticks = 0:0.25:1,
                    ),
                )
            end

            index = findfirst(c -> c.label == curve.label, result.R_T)
            (index === nothing || isempty(result.R_T[index])) && continue
            data = result.R_T[index]
            upper = quantile(data.ratio, 0.99)
            isempty(curve.R_T) || (upper = max(upper, maximum(curve.R_T.ratio)))
            written["figure/R_T_vs_A_H/$(curve.label)"] = FissionTemperatureRatio.save_figure(
                joinpath(directory, "R_T_vs_A_H_segmented_$(token).pdf"),
                FissionTemperatureRatio.plot_ratio(
                    RatioCurve[data],
                    vcat(RatioCurve[curve.R_T], trend_R_T);
                    ylabel = L"R_T = T_L / T_H",
                    reference = 1.0,
                    reference_label = "Equal temperatures",
                    order = order,
                    limits = (first(masses) - 1, last(masses) + 1, 0, 1.15 * upper),
                    annotation = annotation,
                    fit_label = "Segmented fit",
                ),
            )
        end
        return nothing
    end
    return written
end

end
