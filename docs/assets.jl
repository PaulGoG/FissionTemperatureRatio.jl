# Figures and animations that document the method, written into `docs/src/assets/` and shown in
# the README. Regenerate with
#
#     julia docs/assets.jl
#
# Deterministic: the same configuration and input data give the same output, so a regenerated
# asset differs only when the method or the data does.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie
using FissionFragmentsDomain:
    Nuclide, neutron_induced_fission, spontaneous_fission, system_notation
using FissionTemperatureRatio
using LaTeXStrings: @L_str
using Statistics: quantile

const ASSETS = joinpath(@__DIR__, "src", "assets")
const WIDTH = 900
const DATA_DIRECTORY = joinpath(dirname(@__DIR__), "data")
"The multiplicity ratio of the dataset with the most points."
function richest_measurement(configuration)
    A₀ = configuration.system.compound.A
    range = A_H_range(configuration)
    datasets = read_multiplicity_directory(configuration.multiplicity_directory)
    curves = [multiplicity_ratio(data, A₀, range) for data in datasets]
    return curves[argmax(length.(curves))]
end

"""
Animate the model selection across four fissioning systems at once: each panel shows one system's
richest measurement fitted with an increasing number of joined segments, the pivots moving as the
fit gains freedom, and the criterion deciding where each of them stops. The systems do not stop at
the same order, which is the point of showing them together.

Points carry the uncertainty propagated from the multiplicity data, where the archive quotes one.
An order at which no fit stays inside the physical range is shown without one.
"""
function animate_selection(panels; path, max_segments = 6, hold = 8)
    prepared = map(panels) do panel
        fits = map(1:max_segments) do k
            try
                fit_segments(
                    panel.curve.A_H,
                    panel.curve.ratio,
                    panel.curve.σ;
                    min_segments = k,
                    max_segments = k,
                    min_points_per_segment = panel.min_points_per_segment,
                    pinned_value = panel.pinned_value,
                    bounds = (0.0, 1.0),
                )
            catch error
                error isa InsufficientDataError || rethrow()
                nothing
            end
        end
        admissible = findall(!isnothing, fits)
        chosen = admissible[argmin([fits[k].bic for k in admissible])]
        return (; panel.curve, panel.label, fits, chosen)
    end
    # Each panel settles on its own selected order, so the animation ends on four answers.
    frames = vcat(1:max_segments, fill(0, hold))

    figure = Figure(; size = (1200, 900), figure_padding = (10, 20, 10, 10))
    axes = Axis[]
    for (index, prep) in enumerate(prepared)
        row, column = fldmod1(index, 2)
        axis = Axis(
            figure[row, column];
            xlabel = row == 2 ? L"Heavy fragment mass number $A_H$" : "",
            ylabel = column == 1 ? L"r_\nu = \nu_H / (\nu_L + \nu_H)" : "",
            # Every panel keeps its tick labels: the systems span different mass ranges, so
            # hiding them on the top row would imply a shared abscissa that does not exist.
            xticklabelsvisible = true,
            yticklabelsvisible = column == 1,
            xticks = WilkinsonTicks(5),
        )
        xlims!(axis, first(prep.curve.A_H) - 1, last(prep.curve.A_H) + 1)
        ylims!(axis, 0, 1)
        push!(axes, axis)
    end
    colgap!(figure.layout, 10)
    rowgap!(figure.layout, 10)

    record(figure, path, frames; framerate = 2) do frame
        for (axis, prep) in zip(axes, prepared)
            order = frame == 0 ? prep.chosen : min(frame, length(prep.fits))
            fit = prep.fits[order]
            selected = order == prep.chosen
            empty!(axis)
            plural = order == 1 ? "" : "s"

            quoted = findall(!ismissing, prep.curve.σ)
            if !isempty(quoted)
                errorbars!(
                    axis,
                    prep.curve.A_H[quoted],
                    prep.curve.ratio[quoted],
                    Float64[prep.curve.σ[i] for i in quoted];
                    color = (:grey30, 0.5),
                    linewidth = 1.5,
                    whiskerwidth = 0,
                )
            end
            scatter!(
                axis,
                prep.curve.A_H,
                prep.curve.ratio;
                color = (:grey30, 0.6),
                strokecolor = (:grey10, 0.6),
                markersize = 14,
            )

            if fit === nothing
                text!(
                    axis,
                    0.03,
                    0.96;
                    text = "$(prep.label)\n$(order) segment$(plural): none within (0, 1)",
                    space = :relative,
                    align = (:left, :top),
                    fontsize = 21,
                )
                continue
            end
            evaluated = evaluate(fit, first(prep.curve.A_H):last(prep.curve.A_H))
            lines!(
                axis,
                evaluated.A_H,
                evaluated.ratio;
                color = selected ? RGBf(0.0, 0.62, 0.451) : RGBf(0.0, 0.447, 0.698),
                linewidth = 3,
            )
            points = pivots(fit)
            scatter!(
                axis,
                [p[1] for p in points],
                [p[2] for p in points];
                color = :black,
                marker = :diamond,
                markersize = 20,
            )
            text!(
                axis,
                0.03,
                0.96;
                text = "$(prep.label)\n$(order) segment$(plural)" *
                       (selected ? "   ← selected" : ""),
                space = :relative,
                align = (:left, :top),
                fontsize = 21,
            )
        end
        return nothing
    end
    return path
end

# Table 1 of Eur. Phys. J. A 60, 190 (2024), for the comparisons this package can
# make: the data sets and yield distributions it holds. The 239-Pu caption does not name its
# yield distribution; Nishio 1995 is inferred by agreement, and the legend says so.
const PUBLISHED = [
    ("U233_nth", "K. Nishio 1998", "V.M. Surin 1972", 1.1861, 0.0021),
    ("U233_nth", "V.F. Apalin 1965", "V.M. Surin 1972", 1.0346, 0.0043),
    ("U233_nth", "J.S. Fraser 1966", "V.M. Surin 1972", 1.3809, 0.2515),
    ("Cf252_sf", "C. Budtz-jorgensen 1988", "A. Goeoek 2014", 1.0975, 0.0001),
    ("Cf252_sf", "Yu.S. Zamyatnin 1979", "A. Goeoek 2014", 1.1128, 0.0093),
    ("Cf252_sf", "A. Goeoek 2014", "A. Goeoek 2014", 1.1163, 0.0010),
    ("Cf252_sf", "A. Al-adili 2020", "A. Goeoek 2014", 1.0860, 0.0027),
    ("U235_nth", "K. Nishio 1998", "A. Al-adili 2020", 1.1586, 0.0069),
    ("U235_nth", "K. Nishio 1998", "Ch.Straede 1987", 1.1644, 0.0072),
    ("U235_nth", "A.S. Vorobyev 2010", "A. Al-adili 2020", 1.1186, 0.0031),
    ("U235_nth", "A.S. Vorobyev 2010", "Ch.Straede 1987", 1.1221, 0.0031),
    ("Pu239_nth", "K. Nishio 1995", "K. Nishio 1995", 0.9493, 0.0001),
    ("Pu239_nth", "V.F. Apalin 1965", "K. Nishio 1995", 1.0009, 0.0045),
    ("Pu239_nth", "Yu.S. Zamyatnin 1979", "K. Nishio 1995", 1.1588, 0.0146),
    ("Pu239_nth", "J.S. Fraser 1966", "K. Nishio 1995", 1.1072, 0.0149),
]

const SYSTEM_COLOR = Dict(
    "U233_nth" => RGBf(0.0, 0.447, 0.698),
    "Cf252_sf" => RGBf(0.835, 0.369, 0.0),
    "U235_nth" => RGBf(0.0, 0.62, 0.451),
    "Pu239_nth" => RGBf(0.8, 0.475, 0.655),
)
const SYSTEM_MARKER = Dict(
    "U233_nth" => :circle,
    "Cf252_sf" => :rect,
    "U235_nth" => :utriangle,
    "Pu239_nth" => :diamond,
)
const SYSTEM_NOTATION = Dict(
    "U233_nth" => system_notation(neutron_induced_fission(Nuclide(92, 233), 2.53e-8, "nth")),
    "Cf252_sf" => system_notation(spontaneous_fission(Nuclide(98, 252))),
    "U235_nth" => system_notation(neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")),
    "Pu239_nth" =>
        system_notation(neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth")) *
        ", Y(A) inferred",
)

save_asset(name, figure) = save(joinpath(ASSETS, name), figure; px_per_unit = 4)

# Okabe-Ito, the same order the package's own figures use, with seven markers against eight
# colours so that the pairs stay distinct.
const SERIES_COLOURS = [
    RGBf(0.0, 0.447, 0.698),
    RGBf(0.835, 0.369, 0.0),
    RGBf(0.0, 0.62, 0.451),
    RGBf(0.8, 0.475, 0.655),
    RGBf(0.337, 0.706, 0.914),
    RGBf(0.902, 0.624, 0.0),
    RGBf(0.35, 0.35, 0.35),
]
const SERIES_MARKERS = [:circle, :rect, :utriangle, :diamond, :dtriangle, :xcross, :star5]
dataset_colour(i::Integer) = SERIES_COLOURS[mod1(i, length(SERIES_COLOURS))]
# Marker edge: the fill darkened, as the package's own figures do.
darker(colour::RGBf) = RGBf(0.6f0 * colour.r, 0.6f0 * colour.g, 0.6f0 * colour.b)
dataset_marker_symbol(i::Integer) = SERIES_MARKERS[mod1(i, length(SERIES_MARKERS))]

"Published total averages against the ones this package produces, with the residuals beneath."
function figure_published_comparison(results)
    points = NamedTuple[]
    for (system, dataset, yield, published, published_uncertainty) in PUBLISHED
        result = get(results, system, nothing)
        result === nothing && continue
        averages = get(result.total_average_R_T, dataset, nothing)
        averages === nothing && continue
        haskey(averages, yield) || continue
        extracted = averages[yield].value
        extracted_uncertainty = averages[yield].uncertainty
        push!(
            points,
            (; system, published, published_uncertainty, extracted, extracted_uncertainty),
        )
    end

    figure = Figure(; size = (900, 1000))
    axis = Axis(
        figure[1, 1]; ylabel = L"This work, $\langle R_T \rangle$", xticklabelsvisible = false
    )
    residual = Axis(
        figure[2, 1]; xlabel = L"Published $\langle R_T \rangle$", ylabel = "Dev. [%]"
    )
    linkxaxes!(axis, residual)
    rowsize!(figure.layout, 2, Relative(0.26))
    rowgap!(figure.layout, 10)

    low = minimum(p.published for p in points) - 0.03
    high = maximum(p.published for p in points) + 0.03
    lines!(
        axis, [low, high], [low, high]; color = :black, linestyle = :dashdot, linewidth = 1.5
    )
    hlines!(residual, [0.0]; color = :black, linestyle = :dashdot, linewidth = 1.5)
    # The one-per-cent band the comparison is judged against.
    band!(residual, [low, high], [-1.0, -1.0], [1.0, 1.0]; color = (:grey, 0.18))

    for system in ("U233_nth", "Cf252_sf", "U235_nth", "Pu239_nth")
        selected = filter(p -> p.system == system, points)
        isempty(selected) && continue
        x = [p.published for p in selected]
        y = [p.extracted for p in selected]
        errorbars!(
            axis,
            x,
            y,
            [p.extracted_uncertainty for p in selected];
            color = SYSTEM_COLOR[system],
            linewidth = 1.5,
        )
        errorbars!(
            axis,
            x,
            y,
            [p.published_uncertainty for p in selected];
            direction = :x,
            color = SYSTEM_COLOR[system],
            linewidth = 1.5,
        )
        scatter!(
            axis,
            x,
            y;
            color = SYSTEM_COLOR[system],
            strokecolor = darker(SYSTEM_COLOR[system]),
            marker = SYSTEM_MARKER[system],
            label = SYSTEM_NOTATION[system],
        )
        scatter!(
            residual,
            x,
            100 .* (y .- x) ./ x;
            color = SYSTEM_COLOR[system],
            strokecolor = darker(SYSTEM_COLOR[system]),
            marker = SYSTEM_MARKER[system],
        )
    end
    worst = maximum(abs(100 * (p.extracted - p.published) / p.published) for p in points)
    text!(
        axis,
        0.04,
        0.92;
        text = "$(length(points)) comparisons\nlargest deviation $(round(worst; digits = 2)) %",
        space = :relative,
        align = (:left, :top),
        fontsize = 21,
    )
    axislegend(
        axis,
        "Table 1";
        position = :rb,
        framevisible = false,
        labelsize = 22,
        titlesize = 26,
        patchsize = (30, 22),
    )
    # Symmetric about zero, with room above the largest deviation for its marker.
    bound = max(2.0, ceil(worst + 0.4))
    ylims!(residual, -bound, bound)
    # Bounded by the comparisons themselves: Fraser's published uncertainty is a quarter of a
    # unit and would otherwise set the scale for everything else.
    xlims!(axis, low, high)
    ylims!(axis, low, high)
    return figure
end

"The temperature ratio of one system: every measurement, and the systematic trend through them."
function figure_temperature_ratio(result)
    figure = Figure(; size = (900, 650))
    axis = Axis(
        figure[1, 1]; xlabel = L"Heavy fragment mass number $A_H$", ylabel = L"R_T = T_L / T_H"
    )
    hlines!(axis, [1.0]; color = :black, linestyle = :dashdot, linewidth = 1.5)
    for curve in result.R_T
        isempty(curve) && continue
        scatter!(
            axis,
            curve.A_H,
            curve.ratio;
            color = (:grey45, 0.55),
            strokecolor = (:grey25, 0.55),
            markersize = 14,
        )
    end
    trend = systematic_trend(result).R_T
    band!(
        axis,
        trend.A_H,
        max.(trend.ratio .- trend.σ, 0.0),
        trend.ratio .+ trend.σ;
        color = (RGBf(0.835, 0.369, 0.0), 0.25),
    )
    lines!(axis, trend.A_H, trend.ratio; color = RGBf(0.835, 0.369, 0.0), linewidth = 3)
    # Every measurement is drawn; the trend is fitted to those the configuration pools.
    excluded = result.configuration.excluded_datasets
    pooled = count(data -> !haskey(excluded, data.label), result.datasets)
    through = pooled == length(result.datasets) ? "" : " through $(pooled)"
    text!(
        axis,
        0.97,
        0.93;
        text = "$(system_notation(result.configuration.system))\n\
                $(length(result.datasets)) measurements\nSystematic trend$(through)",
        space = :relative,
        align = (:right, :top),
        fontsize = 21,
    )
    ylims!(axis, 0, 3)
    return figure
end

"""
From the measured sawtooth to the temperature ratio, on a shared abscissa, for several
measurements at once — and the curve combined from all of them.

Several rather than one, because the step that needs showing is how measurements of the same
quantity are used together: each is carried through separately, and only the combined curve of the
middle panel merges them.
"""
function figure_method_chain(result, count = 5)
    A₀ = result.configuration.system.compound.A
    order = sortperm(length.(result.r_ν); rev = true)
    chosen = [i for i in order if !isempty(result.r_ν[i])][1:min(count, length(order))]

    figure = Figure(; size = (900, 1300))
    axes = [
        Axis(figure[row, 1]; ylabel = ylabel, xticklabelsvisible = row == 3) for
        (row, ylabel) in enumerate((L"\nu", L"r_\nu = \nu_H/(\nu_L + \nu_H)", L"R_T = T_L/T_H"))
    ]
    axes[3].xlabel = L"Fragment mass number $A$"
    linkxaxes!(axes...)
    # Wide enough that the lowest tick label of one panel clears the highest of the next: each
    # panel is bounded at a round value, so both labels sit on the frame and a narrower gap runs
    # them together.
    rowgap!(figure.layout, 28)

    # The symmetric split separates the light wing from the heavy one, which is what a label on
    # each wing used to say — and it says it without sitting on top of the data.
    for axis in axes
        vlines!(axis, [A₀ / 2]; color = (:grey, 0.55), linestyle = :dot, linewidth = 1.5)
    end
    hlines!(axes[2], [0.5]; color = :black, linestyle = :dashdot, linewidth = 1.5)
    hlines!(axes[3], [1.0]; color = :black, linestyle = :dashdot, linewidth = 1.5)

    for (position, index) in enumerate(chosen)
        data = result.datasets[index]
        colour = dataset_colour(position)
        marker = dataset_marker_symbol(position)
        scatter!(
            axes[1],
            data.A,
            data.ν;
            color = colour,
            strokecolor = darker(colour),
            marker = marker,
            markersize = 14,
            label = data.label,
        )
        for (axis, curve) in ((axes[2], result.r_ν[index]), (axes[3], result.R_T[index]))
            scatter!(
                axis,
                curve.A_H,
                curve.ratio;
                color = colour,
                strokecolor = darker(colour),
                marker = marker,
                markersize = 14,
            )
        end
    end

    # The combined curve, which is the only place the measurements are merged.
    consensus = result.consensus_r_ν
    lines!(
        axes[2],
        consensus.A_H,
        consensus.ratio;
        color = :black,
        linestyle = :dash,
        linewidth = 3,
        label = "combined",
    )
    trend = systematic_trend(result).R_T
    lines!(axes[3], trend.A_H, trend.ratio; color = :black, linestyle = :dash, linewidth = 3)

    Legend(
        figure[0, 1],
        axes[1];
        orientation = :horizontal,
        nbanks = 3,
        framevisible = false,
        labelsize = 22,
        patchsize = (30, 22),
        tellheight = true,
        tellwidth = false,
        padding = (0, 0, 0, 0),
    )
    # Set after the legend, because adding a row renumbers the gaps: indexing one of them by hand
    # tightened the wrong one and ran two panels' tick labels together.
    rowgap!(figure.layout, 28)
    # The far-asymmetric tail of some measurements reaches tens of neutrons per fragment; scaled
    # to those, the sawtooth that carries the physics collapses to a line. Bounded by the bulk, as
    # the package's own multiplicity figure is.
    # The mass range the analysis actually uses, symmetric about the split.
    masses = A_H_range(result.configuration)
    low, high = A₀ - last(masses) - 2, last(masses) + 2
    xlims!(axes[1], low, high)
    bulk = reduce(
        vcat,
        (result.datasets[i].ν[low .≤ result.datasets[i].A .≤ high] for i in chosen);
        init = Float64[],
    )
    ylims!(axes[1], 0, 1.1 * maximum(bulk))
    ylims!(axes[2], 0, 1)
    ylims!(axes[3], 0, 3)
    text!(
        axes[1],
        0.02,
        0.94;
        text = system_notation(result.configuration.system),
        space = :relative,
        align = (:left, :top),
        fontsize = 21,
        color = :grey40,
    )
    return figure
end

"How much the level density model moves the answer."
function figure_level_density_models(bsfg, gc)
    figure = Figure(; size = (900, 650))
    axis = Axis(
        figure[1, 1]; xlabel = L"Heavy fragment mass number $A_H$", ylabel = L"R_T = T_L / T_H"
    )
    hlines!(axis, [1.0]; color = :black, linestyle = :dashdot, linewidth = 1.5)
    for (result, name, colour, style) in (
        (bsfg, "Back-shifted Fermi gas", RGBf(0.0, 0.447, 0.698), :solid),
        (gc, "Gilbert-Cameron", RGBf(0.835, 0.369, 0.0), :dash),
    )
        trend = systematic_trend(result).R_T
        lines!(
            axis,
            trend.A_H,
            trend.ratio;
            color = colour,
            linewidth = 3,
            linestyle = style,
            label = name,
        )
    end
    axislegend(axis; position = :rt, framevisible = false, labelsize = 22, patchsize = (30, 22))
    text!(
        axis,
        0.03,
        0.05;
        text = "$(system_notation(bsfg.configuration.system))\n\
                systematic trend under each level density model",
        space = :relative,
        align = (:left, :bottom),
        fontsize = 21,
        color = :grey40,
    )
    return figure
end

# A system under its shipped settings, averaged over every yield distribution held rather than the
# primary experiment's alone, which stays the coverage reference: Table 1 names distributions the
# shipped configurations do not take.
function every_distribution(system)
    source = joinpath(dirname(@__DIR__), "config", "$(system).toml")
    text = replace(
        read(source, String),
        r"^(mass_yield_file = .*)$"m =>
            SubstitutionString("\\1\nsubdirectory = \"$(system)/Y_vs_A\""),
    )
    return mktempdir() do directory
        path = joinpath(directory, "$(system).toml")
        write(path, text)
        return run_pipeline(load_configuration(path; data_directory = DATA_DIRECTORY))
    end
end

function main()
    mkpath(ASSETS)
    systems = ("U233_nth", "U235_nth", "Pu239_nth", "Cf252_sf")
    configurations = Dict(
        system => load_configuration(
            joinpath(dirname(@__DIR__), "config", "$(system).toml");
            data_directory = DATA_DIRECTORY,
        ) for system in systems
    )

    panels = map(systems) do system
        configuration = configurations[system]
        curve = richest_measurement(configuration)
        settings = configuration.segments
        # As the pipeline does: the pin is an identity at the symmetric split and nowhere else.
        pinned =
            if settings.pin_symmetric_split &&
                2 * first(curve.A_H) == configuration.system.compound.A
                0.5
            else
                nothing
            end
        return (;
            curve,
            label = "$(system_notation(configuration.system)) · $(curve.label)",
            min_points_per_segment = settings.min_points_per_segment,
            pinned_value = pinned,
        )
    end

    with_theme(publication_theme()) do
        path = animate_selection(
            panels;
            path = joinpath(ASSETS, "segment_selection.gif"),
            max_segments = maximum(configurations[s].segments.max_segments for s in systems),
        )
        @info "written" path filesize = filesize(path)
    end

    # The static figures need the pipeline run rather than one fitted curve. Output is not written
    # again here; these read the results in memory.
    results = Dict{String,ExtractionResult}()
    for system in ("U233_nth", "U235_nth", "Cf252_sf")
        results[system] = run_pipeline(configurations[system])
    end
    # The same system under the other level density model, for the comparison of the two.
    source = joinpath(dirname(@__DIR__), "config", "U233_nth.toml")
    gilbert_cameron = mktempdir() do directory
        path = joinpath(directory, "gilbert_cameron.toml")
        write(path, replace(read(source, String), "model = \"BSFG\"" => "model = \"GC\""))
        return run_pipeline(load_configuration(path; data_directory = DATA_DIRECTORY))
    end

    with_theme(publication_theme()) do
        for (name, figure) in (
            (
                "published_comparison.png",
                figure_published_comparison(
                    Dict(
                        s => every_distribution(s) for
                        s in ("U233_nth", "U235_nth", "Cf252_sf", "Pu239_nth")
                    ),
                ),
            ),
            ("temperature_ratio.png", figure_temperature_ratio(results["Cf252_sf"])),
            ("method_chain.png", figure_method_chain(results["Cf252_sf"])),
            (
                "level_density_models.png",
                figure_level_density_models(results["U233_nth"], gilbert_cameron),
            ),
        )
            save_asset(name, figure)
            @info "written" name filesize = filesize(joinpath(ASSETS, name))
        end
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
