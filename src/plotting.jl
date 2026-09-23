# Publication figures.
#
# Only the interface is declared here. The implementations are in the CairoMakie extension, so
# that a package whose output is tabulated data does not load a plotting stack on `using`. Loading
# CairoMakie supplies them:
#
#     using FissionTemperatureRatio, CairoMakie
#
# Without it these functions exist but have no methods; the tables and the manifest a run writes
# do not need them.

"""
    publication_theme() -> Theme

The figure style used throughout, the standard layout: Computer Modern faces, a 900 × 600 canvas
per panel, 26 pt type, 3 pt data lines, 14 pt markers with a darker edge, boxed axes with inward
ticks, no minor ticks, faint dashed gridlines.

Requires CairoMakie to be loaded.
"""
function publication_theme end

"""
    plot_multiplicities(datasets; A_0) -> Figure

Prompt neutron multiplicity `ν(A)` of every dataset, the sawtooth that carries the signature of
the excitation energy partition.

Requires CairoMakie to be loaded.
"""
function plot_multiplicities end

"""
    plot_ratio(curves, fitted; ylabel, reference, reference_label, annotation, annotation_corner, fit_label,
               yticks) -> Figure

A ratio against the heavy-fragment mass number: the values extracted from each dataset with
their uncertainties, and optionally a fitted curve with its uncertainty band.

`reference` draws a horizontal guide line at a value the ratio takes by construction — one half
for the multiplicity ratio, unity for the temperature ratio — labelled so that the reader need
not infer what it marks.

`annotation` puts the quantitative takeaway inside the axes, at the corner `annotation_corner`
names: `:rt`, the default, which a temperature ratio leaves free, or `:rb`, which a multiplicity
ratio does. Give it one
string per line, or a single string for one line; each line is typeset on its own, so a
`LaTeXString` line renders as mathematics rather than as its own source.

`fit_label` names a dataset's own fitted curve in the legend; the systematic-trend curve is always
labelled, and is drawn as a thin dashed guide when it accompanies another fitted curve.

Requires CairoMakie to be loaded.
"""
function plot_ratio end

"""
    save_figure(path, figure) -> String

Write a figure to `path` as vector PDF, creating the directory if needed and never overwriting an
existing file. Returns the path actually written.

Requires CairoMakie to be loaded.
"""
function save_figure end

"""
    write_figures(result, directory) -> Dict{String,String}

Write the figures of a completed run into `directory`, returning the paths written, keyed by
content. Besides the three overview figures, every per-dataset segmented curve gets its own
multiplicity-ratio and temperature-ratio figure, with the systematic trend as a guide. The
directory names the run; the file names name the quantity and the abscissa.

Requires CairoMakie to be loaded. `scripts/run.jl` calls it after [`write_results`](@ref), into
`plots/<system>/<run>/`.
"""
function write_figures end
