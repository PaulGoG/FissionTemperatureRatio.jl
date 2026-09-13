# Figures and animations that document the method, written into `docs/src/assets/` and shown in
# the README. Regenerate with
#
#     julia --project=docs docs/assets.jl [config]
#
# Deterministic: the same configuration and input data give the same output, so a regenerated
# asset differs only when the method or the data does.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie
using FissionTemperatureRatio
using LaTeXStrings: @L_str

const ASSETS = joinpath(@__DIR__, "src", "assets")
const CONFIG =
    length(ARGS) ≥ 1 ? ARGS[1] : joinpath(dirname(@__DIR__), "config", "Cf252_0f.toml")

"The multiplicity ratio of the data set with the most points, and the level density ratio."
function reference_case(configuration)
    prescription = build_prescription(configuration.level_density)
    charges = read_charge_distribution(
        configuration.fragmentation.charge_distribution_file,
        configuration.fragmentation.fallback_charge_polarization,
        configuration.fragmentation.fallback_rms,
    )
    A₀, Z₀ = configuration.system.A₀, configuration.system.Z₀
    range = A_H_range(configuration)
    domain = fragmentation_domain(
        A₀, Z₀, range, configuration.fragmentation.charges_per_mass, charges
    )
    R_a = level_density_ratio(
        configuration.level_density.ratio_averaging, prescription, A₀, Z₀, domain
    )
    data = read_multiplicity_directory(configuration.multiplicity_directory)
    curves = [multiplicity_ratio(set, A₀, range) for set in data]
    return (curves[argmax(length.(curves))], R_a)
end

"""
Animate the model selection: the same data fitted with an increasing number of segments, the
pivots moving as the fit gains freedom, and the criterion deciding where to stop.
"""
function animate_selection(curve, settings; path, caption, max_segments = 6, hold = 8)
    orders = vcat(1:max_segments, fill(max_segments, 0))
    fits = [
        fit_segments(
            curve.A_H,
            curve.value,
            curve.σ;
            min_segments = k,
            max_segments = k,
            min_points_per_segment = settings.min_points_per_segment,
            pinned_value = settings.pinned_value,
            bounds = (0.0, 1.0),
        ) for k in orders
    ]
    chosen = argmin(fit.bic for fit in fits)
    # The chosen order is held at the end, so the animation settles on the answer.
    frames = vcat(1:length(fits), fill(chosen, hold))

    # Extra padding on the right: the last tick label sits on the frame and is otherwise clipped.
    figure = Figure(; size = (470, 300), figure_padding = (6, 14, 6, 8))
    axis = Axis(
        figure[1, 1];
        xlabel = L"Heavy fragment mass number $A_H$",
        ylabel = L"r_\nu = \nu_H / (\nu_L + \nu_H)",
        xticks = WilkinsonTicks(6),
    )
    xlims!(axis, first(curve.A_H) - 1, last(curve.A_H) + 1)
    ylims!(axis, 0, 1)

    record(figure, path, frames; framerate = 2) do index
        fit = fits[index]
        empty!(axis)
        scatter!(
            axis,
            curve.A_H,
            curve.value;
            color = (:grey30, 0.55),
            markersize = 5,
            label = curve.label,
        )
        evaluated = evaluate(fit, first(curve.A_H):last(curve.A_H))
        selected = index == chosen
        lines!(
            axis,
            evaluated.A_H,
            evaluated.value;
            color = selected ? RGBf(0.0, 0.62, 0.451) : RGBf(0.0, 0.447, 0.698),
            linewidth = 2,
        )
        points = pivots(fit)
        scatter!(
            axis,
            [p[1] for p in points],
            [p[2] for p in points];
            color = :black,
            marker = :diamond,
            markersize = 11,
        )
        text!(
            axis,
            0.03,
            0.95;
            text = "$(segments(fit)) segment$(segments(fit) == 1 ? "" : "s")\nBIC $(round(Int, fit.bic))" *
                   (selected ? "   ← selected" : ""),
            space = :relative,
            align = (:left, :top),
            fontsize = 13,
        )
        text!(
            axis,
            0.97,
            0.05;
            text = caption,
            space = :relative,
            align = (:right, :bottom),
            fontsize = 11,
            color = :grey40,
        )
        return nothing
    end
    return path
end

function main()
    mkpath(ASSETS)
    # The active project here is `docs/`, so DrWatson's default data directory would resolve
    # under it; the data lives beside the package.
    configuration = load_configuration(
        CONFIG; data_directory = joinpath(dirname(@__DIR__), "data")
    )
    curve, _ = reference_case(configuration)
    settings = configuration.segments
    pinned = settings.pin_symmetric_split && has_symmetric_split(configuration) ? 0.5 : nothing

    with_theme(publication_theme(); fontsize = 13) do
        path = animate_selection(
            curve,
            (min_points_per_segment = settings.min_points_per_segment, pinned_value = pinned);
            path = joinpath(ASSETS, "segment_selection.gif"),
            caption = "$(configuration.system.label)   ·   $(curve.label)",
            max_segments = settings.max_segments,
        )
        @info "written" path filesize = filesize(path)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
