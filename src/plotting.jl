# Publication figures. Designed at the printed width of a single journal column, so that they
# enter a manuscript at native size without rescaling.

"""
    SINGLE_COLUMN_WIDTH

Width of a single journal column in points, 86 mm. Figures are laid out at their final printed
size; anything that needs rescaling on inclusion was exported at the wrong size.
"""
const SINGLE_COLUMN_WIDTH = 86 / 25.4 * 72

"""
    publication_theme() -> Theme

The figure style used throughout: Computer Modern faces, boxed axes with inward ticks, no minor
ticks, faint dashed gridlines, and type sizes that remain legible at a single column width.
"""
function publication_theme()
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

# One colour per data set, consistent across every figure of a run. Okabe-Ito, which stays
# distinguishable in grayscale and for the common colour vision deficiencies.
const DATA_SET_COLORS = [
    RGBf(0.0, 0.447, 0.698),
    RGBf(0.835, 0.369, 0.0),
    RGBf(0.0, 0.620, 0.451),
    RGBf(0.800, 0.475, 0.655),
    RGBf(0.337, 0.706, 0.914),
    RGBf(0.941, 0.894, 0.259),
    RGBf(0.902, 0.624, 0.0),
    RGBf(0.35, 0.35, 0.35),
]

const DATA_SET_MARKERS = [
    :circle, :rect, :utriangle, :diamond, :dtriangle, :cross, :xcross, :star5
]

data_set_color(index::Integer) = DATA_SET_COLORS[mod1(index, length(DATA_SET_COLORS))]
data_set_marker(index::Integer) = DATA_SET_MARKERS[mod1(index, length(DATA_SET_MARKERS))]

"""
    plot_multiplicities(data_sets; A₀) -> Figure

Prompt neutron multiplicity `ν(A)` of every data set, the sawtooth that carries the signature of
the excitation energy partition.
"""
function plot_multiplicities(data_sets::Vector{MultiplicityData}; A₀::Integer)
    figure = Figure(; size = (SINGLE_COLUMN_WIDTH, 0.78 * SINGLE_COLUMN_WIDTH))
    axis = Axis(
        figure[1, 1];
        xlabel = L"Fragment mass number $A$",
        ylabel = L"Prompt neutron multiplicity $\nu$",
    )

    for (index, data) in enumerate(data_sets)
        color = data_set_color(index)
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
            marker = data_set_marker(index),
            label = data.label,
        )
    end

    _legend_above(figure, axis, length(data_sets))
    return figure
end

"""
    plot_ratio(curves, fitted; ylabel, reference, reference_label, annotation) -> Figure

A ratio against the heavy-fragment mass number: the values extracted from each data set with
their uncertainties, and optionally a fitted curve with its uncertainty band.

`reference` draws a horizontal guide line at a value the ratio takes by construction — one half
for the multiplicity ratio, unity for the temperature ratio — labelled so that the reader need
not infer what it marks.
"""
function plot_ratio(
    curves::Vector{RatioCurve},
    fitted::Union{RatioCurve,Nothing} = nothing;
    ylabel,
    reference::Union{Real,Nothing} = nothing,
    reference_label::AbstractString = "",
    annotation::AbstractString = "",
)
    figure = Figure(; size = (SINGLE_COLUMN_WIDTH, 0.78 * SINGLE_COLUMN_WIDTH))
    axis = Axis(figure[1, 1]; xlabel = L"Heavy fragment mass number $A_H$", ylabel = ylabel)

    if reference !== nothing
        hlines!(axis, [reference]; color = :black, linestyle = :dashdot, linewidth = 0.7)
    end

    for (index, curve) in enumerate(curves)
        isempty(curve) && continue
        color = data_set_color(index)
        if any(>(0), curve.σ)
            errorbars!(
                axis,
                curve.A_H,
                curve.value,
                curve.σ;
                color = color,
                linewidth = 0.6,
                whiskerwidth = 3,
            )
        end
        scatter!(
            axis,
            curve.A_H,
            curve.value;
            color = color,
            marker = data_set_marker(index),
            label = curve.label,
        )
    end

    if fitted !== nothing && !isempty(fitted)
        band!(
            axis,
            fitted.A_H,
            fitted.value .- fitted.σ,
            fitted.value .+ fitted.σ;
            color = (:black, 0.15),
        )
        lines!(
            axis,
            fitted.A_H,
            fitted.value;
            color = :black,
            linewidth = 1.2,
            label = fitted.label,
        )
    end

    if !isempty(reference_label) && reference !== nothing
        # Both coordinates are in data space: the label marks a particular value of the ratio, so
        # its ordinate is that value and not a fraction of the axis height.
        left = minimum(
            minimum(curve.A_H) for curve in curves if !isempty(curve);
            init = fitted === nothing ? 0 : minimum(fitted.A_H),
        )
        text!(
            axis,
            left,
            reference;
            text = reference_label,
            align = (:left, :bottom),
            fontsize = 6,
            offset = (0, 2),
        )
    end
    if !isempty(annotation)
        text!(
            axis,
            0.97,
            0.05;
            text = annotation,
            space = :relative,
            align = (:right, :bottom),
            fontsize = 6,
        )
    end

    _legend_above(figure, axis, count(!isempty, curves) + (fitted === nothing ? 0 : 1))
    return figure
end

# Legends sit above the axes, horizontally, so that they cannot collide with the data however the
# points and their error bars happen to fall.
function _legend_above(figure::Figure, axis::Axis, entries::Integer)
    # At a single column width three entries fit on one line; more must be banked, or the last
    # label is silently clipped.
    Legend(
        figure[0, 1],
        axis;
        orientation = :horizontal,
        nbanks = entries ≤ 3 ? 1 : cld(entries, 3),
        framevisible = false,
        labelsize = 6,
        padding = (0, 0, 0, 0),
        tellheight = true,
        tellwidth = false,
    )
    rowgap!(figure.layout, 2)
    return figure
end

"""
    save_figure(path, figure)

Write a figure to `path` as vector PDF, creating the directory if needed and never overwriting an
existing file. Returns the path actually written.
"""
function save_figure(path::AbstractString, figure::Figure)
    mkpath(dirname(path))
    target = _unused_path(path)
    save(target, figure; pt_per_unit = 1)
    return target
end
