# Publication figures.
#
# Only the interface is declared here. The implementations are in the CairoMakie extension, so
# that a package whose output is tabulated data does not load a plotting stack on `using`. Loading
# CairoMakie supplies them:
#
#     using FissionTemperatureRatio, CairoMakie
#
# Without it these functions exist but have no methods, and `run_pipeline` writes its tables and
# metadata while reporting that figures were skipped.

"""
    publication_theme() -> Theme

The figure style used throughout: Computer Modern faces, boxed axes with inward ticks, no minor
ticks, faint dashed gridlines, and type sizes that remain legible at a single column width.

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
    plot_ratio(curves, fitted; ylabel, reference, reference_label) -> Figure

A ratio against the heavy-fragment mass number: the values extracted from each dataset with
their uncertainties, and optionally a fitted curve with its uncertainty band.

`reference` draws a horizontal guide line at a value the ratio takes by construction — one half
for the multiplicity ratio, unity for the temperature ratio — labelled so that the reader need
not infer what it marks.

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
    write_figures(result, directory, identifier) -> Dict{String,String}

Write the figures of a completed run into `directory`, returning the paths written, keyed by
content.

Requires CairoMakie to be loaded; [`write_results`](@ref) calls this only when it is, and records
that figures were skipped when it is not.
"""
function write_figures end

# Whether the plotting extension is available. `write_results` consults this rather than calling
# into the extension blindly, so that a run without CairoMakie still produces its tables.
function _plotting_extension()
    return Base.get_extension(@__MODULE__, :FissionTemperatureRatioCairoMakieExt)
end
