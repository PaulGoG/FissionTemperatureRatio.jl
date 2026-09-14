# Publication figures, laid out at the printed width of a single journal column so that they enter
# a manuscript at native size without rescaling.
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

const SINGLE_COLUMN_WIDTH = 86 / 25.4 * 72

function FissionTemperatureRatio.publication_theme()
    return Theme(;
        fonts = (; regular = texfont(:text), bold = texfont(:bold), italic = texfont(:italic)),
        fontsize = 8,
        figure_padding = 4,
        Axis = (
            xgridstyle = :dash,
            ygridstyle = :dash,
            xgridcolor = (:grey, 0.12),
            ygridcolor = (:grey, 0.12),
            xminorticksvisible = false,
            yminorticksvisible = false,
            xtickalign = 1,
            ytickalign = 1,
            spinewidth = 0.8,
            xtickwidth = 0.8,
            ytickwidth = 0.8,
            xlabelpadding = 2,
            ylabelpadding = 2,
        ),
        Legend = (
            framevisible = false,
            padding = (2, 2, 2, 2),
            rowgap = 0,
            colgap = 6,
            patchsize = (10, 6),
        ),
        Scatter = (markersize = 4, strokewidth = 0),
        Lines = (linewidth = 1,),
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

# Two legend entries per row: at a single column width three columns of author-and-year labels
# overflow the figure and the rightmost is silently clipped.
const LEGEND_COLUMNS = 2
legend_rows(entries::Integer) = max(1, cld(entries, LEGEND_COLUMNS))

# The axes keep a constant height; the figure grows to make room for the legend above them.
function _figure_size(entries::Integer)
    return (SINGLE_COLUMN_WIDTH, 0.78 * SINGLE_COLUMN_WIDTH + 7.5 * legend_rows(entries))
end

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
        if any(>(0), data.σν)
            errorbars!(
                axis, data.A, data.ν, data.σν; color = color, linewidth = 0.6, whiskerwidth = 3
            )
        end
        scatter!(
            axis,
            data.A,
            data.ν;
            color = color,
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
)
    entries =
        count(!isempty, curves) +
        (reference === nothing ? 0 : 1) +
        count(curve -> !isempty(curve) && curve.label == TREND_LABEL, fitted)
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
            linestyle = :dashdot,
            linewidth = 0.7,
            label = isempty(reference_label) ? nothing : reference_label,
        )
    end

    for (position, curve) in enumerate(curves)
        isempty(curve) && continue
        index = _style_index(curve.label, order, position)
        color = dataset_color(index)
        if any(>(0), curve.σ)
            errorbars!(
                axis,
                curve.A_H,
                curve.ratio,
                curve.σ;
                color = color,
                linewidth = 0.6,
                whiskerwidth = 3,
            )
        end
        scatter!(
            axis,
            curve.A_H,
            curve.ratio;
            color = color,
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
        # Both ratios drawn here are non-negative by construction, so the band is clipped at zero
        # rather than drawn into a region the quantity cannot occupy. The symmetric interval is a
        # Gaussian approximation; where it reaches below zero it is the approximation failing, not
        # the quantity.
        band!(
            axis,
            curve.A_H,
            max.(curve.ratio .- curve.σ, 0.0),
            curve.ratio .+ curve.σ;
            color = (color, 0.15),
        )
        lines!(
            axis,
            curve.A_H,
            curve.ratio;
            color = color,
            linewidth = 1.2,
            linestyle = trend ? :dash : :solid,
            label = trend ? curve.label : nothing,
        )
    end

    # The takeaway belongs where the reader is looking, inside the axes.
    _annotate!(axis, annotation)

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
function _annotate!(axis::Axis, annotation; fontsize = 7, leading = 1.3)
    lines =
        annotation isa AbstractString ? (isempty(annotation) ? () : (annotation,)) : annotation
    for (position, line) in enumerate(lines)
        text!(
            axis,
            0.98,
            0.96;
            text = line,
            space = :relative,
            align = (:right, :top),
            offset = (0, -(position - 1) * leading * fontsize),
            fontsize = fontsize,
        )
    end
    return axis
end

# Legends sit above the axes, horizontally, so that they cannot collide with the data however the
# points and their error bars happen to fall.
function _legend_above(figure::Figure, axis::Axis, entries::Integer)
    # Banked so that no row exceeds LEGEND_COLUMNS entries; wider rows overflow the figure and the
    # rightmost label is clipped without any error.
    Legend(
        figure[0, 1],
        axis;
        orientation = :horizontal,
        nbanks = legend_rows(entries),
        framevisible = false,
        labelsize = 6,
        padding = (0, 0, 0, 0),
        tellheight = true,
        tellwidth = false,
    )
    rowgap!(figure.layout, 2)
    return figure
end

function FissionTemperatureRatio.save_figure(path::AbstractString, figure::Figure)
    mkpath(dirname(path))
    target = FissionTemperatureRatio._unused_path(path)
    save(target, figure; pt_per_unit = 1)
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
        value, uncertainty = averages[name]
        return [
            L"Trend $\langle R_T \rangle = %$(_measured(value, uncertainty))$",
            L"over $Y(A)$ of %$(name)",
        ]
    end
    value, uncertainty = result.range_mean_R_T[TREND_LABEL]
    return [L"Trend range mean $= %$(_measured(value, uncertainty))$"]
end

# The run's figures. Held here, rather than in the pipeline, so that the pipeline carries no
# reference to a plotting type; `write_results` calls it through `Base.get_extension`.
function FissionTemperatureRatio.write_figures(
    result::ExtractionResult, directory::AbstractString, identifier::AbstractString
)
    written = Dict{String,String}()
    configuration = result.configuration
    order = [data.label for data in result.datasets]
    masses = FissionTemperatureRatio.A_H_range(configuration)
    with_theme(FissionTemperatureRatio.publication_theme()) do
        written["figure/nu_vs_A"] = FissionTemperatureRatio.save_figure(
            joinpath(directory, "nu_vs_A_$(identifier).pdf"),
            FissionTemperatureRatio.plot_multiplicities(
                result.datasets; A_0 = configuration.system.A₀, order = order
            ),
        )
        written["figure/r_nu_vs_A_H"] = FissionTemperatureRatio.save_figure(
            joinpath(directory, "r_nu_vs_A_H_$(identifier).pdf"),
            FissionTemperatureRatio.plot_ratio(
                filter(!isempty, result.r_ν),
                [c.r_ν for c in result.segmented_curves];
                ylabel = L"r_\nu = \nu_H / (\nu_L + \nu_H)",
                reference = 0.5,
                reference_label = "Equal sharing",
                order = order,
                limits = (first(masses) - 1, last(masses) + 1, 0, 1),
            ),
        )
        return written["figure/R_T_vs_A_H"] = FissionTemperatureRatio.save_figure(
            joinpath(directory, "R_T_vs_A_H_$(identifier).pdf"),
            FissionTemperatureRatio.plot_ratio(
                filter(!isempty, result.R_T),
                [c.R_T for c in result.segmented_curves];
                ylabel = L"R_T = T_L / T_H",
                reference = 1.0,
                reference_label = "Equal temperatures",
                order = order,
                limits = (first(masses) - 1, last(masses) + 1, nothing, nothing),
                annotation = _trend_annotation(result),
            ),
        )
    end
    return written
end

end
