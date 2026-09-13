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
const SINGLE_COLUMN = 86 / 25.4 * 72
const DATA_DIRECTORY = joinpath(dirname(@__DIR__), "data")
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

# Table 1 and Table 2 of Eur. Phys. J. A 60, 190 (2024), for the comparisons this package can
# make: the data sets and yield distributions it holds. 239-Pu is excluded because the table
# averages over a yield distribution its caption does not name.
const PUBLISHED = [
    ("U233_nf", "K. Nishio 1998", "V.M. Surin 1972", 1.1861, 0.0021),
    ("U233_nf", "V.F. Apalin 1965", "V.M. Surin 1972", 1.0346, 0.0043),
    ("U233_nf", "J.S. Fraser 1966", "V.M. Surin 1972", 1.3809, 0.2515),
    ("Cf252_0f", "C. Budtz-jorgensen 1988", "A. Goeoek 2014", 1.0975, 0.0001),
    ("Cf252_0f", "Yu.S. Zamyatnin 1979", "A. Goeoek 2014", 1.1128, 0.0093),
    ("Cf252_0f", "A. Goeoek 2014", "A. Goeoek 2014", 1.1163, 0.0010),
    ("Cf252_0f", "A. Al-adili 2020", "A. Goeoek 2014", 1.0860, 0.0027),
    ("U235_nf", "K. Nishio 1998", "A. Al-adili 2020", 1.1586, 0.0069),
    ("U235_nf", "K. Nishio 1998", "Ch.Straede 1987", 1.1644, 0.0072),
    ("U235_nf", "A.S. Vorobyev 2010", "A. Al-adili 2020", 1.1186, 0.0031),
    ("U235_nf", "A.S. Vorobyev 2010", "Ch.Straede 1987", 1.1221, 0.0031),
]

const SYSTEM_COLOR = Dict(
    "U233_nf" => RGBf(0.0, 0.447, 0.698),
    "Cf252_0f" => RGBf(0.835, 0.369, 0.0),
    "U235_nf" => RGBf(0.0, 0.62, 0.451),
)
const SYSTEM_MARKER = Dict("U233_nf" => :circle, "Cf252_0f" => :rect, "U235_nf" => :utriangle)

save_asset(name, figure) = save(joinpath(ASSETS, name), figure; px_per_unit = 4)

"Published total averages against the ones this package produces, with the residuals beneath."
function figure_published_comparison(results)
    points = NamedTuple[]
    for (system, set, yield, value, uncertainty) in PUBLISHED
        result = get(results, system, nothing)
        result === nothing && continue
        averages = get(result.total_average_R_T, set, nothing)
        averages === nothing && continue
        haskey(averages, yield) || continue
        ours, σ = averages[yield]
        push!(points, (; system, value, uncertainty, ours, σ))
    end

    figure = Figure(; size = (SINGLE_COLUMN, 1.05 * SINGLE_COLUMN))
    axis = Axis(
        figure[1, 1]; ylabel = L"This work, $\langle R_T \rangle$", xticklabelsvisible = false
    )
    residual = Axis(
        figure[2, 1]; xlabel = L"Published $\langle R_T \rangle$", ylabel = "Dev. [%]"
    )
    linkxaxes!(axis, residual)
    rowsize!(figure.layout, 2, Relative(0.26))
    rowgap!(figure.layout, 6)

    low = minimum(p.value for p in points) - 0.03
    high = maximum(p.value for p in points) + 0.03
    lines!(
        axis, [low, high], [low, high]; color = :black, linestyle = :dashdot, linewidth = 0.7
    )
    hlines!(residual, [0.0]; color = :black, linestyle = :dashdot, linewidth = 0.7)
    # The one-per-cent band the comparison is judged against.
    band!(residual, [low, high], [-1.0, -1.0], [1.0, 1.0]; color = (:grey, 0.18))

    for system in ("U233_nf", "Cf252_0f", "U235_nf")
        selected = filter(p -> p.system == system, points)
        isempty(selected) && continue
        x = [p.value for p in selected]
        y = [p.ours for p in selected]
        errorbars!(
            axis, x, y, [p.σ for p in selected]; color = SYSTEM_COLOR[system], linewidth = 0.6
        )
        errorbars!(
            axis,
            x,
            y,
            [p.uncertainty for p in selected];
            direction = :x,
            color = SYSTEM_COLOR[system],
            linewidth = 0.6,
        )
        scatter!(
            axis,
            x,
            y;
            color = SYSTEM_COLOR[system],
            marker = SYSTEM_MARKER[system],
            label = replace(system, "_nf" => "(n,f)", "_0f" => "(sf)"),
        )
        scatter!(
            residual,
            x,
            100 .* (y .- x) ./ x;
            color = SYSTEM_COLOR[system],
            marker = SYSTEM_MARKER[system],
        )
    end
    worst = maximum(abs(100 * (p.ours - p.value) / p.value) for p in points)
    text!(
        axis,
        0.04,
        0.92;
        text = "$(length(points)) comparisons\nall within $(round(worst; digits = 1)) %",
        space = :relative,
        align = (:left, :top),
        fontsize = 7,
    )
    axislegend(
        axis,
        "Tables 1 and 2";
        position = :rb,
        framevisible = false,
        labelsize = 6.5,
        titlesize = 7,
        patchsize = (8, 6),
    )
    ylims!(residual, -2.0, 2.0)
    # Bounded by the comparisons themselves: Fraser's published uncertainty is a quarter of a
    # unit and would otherwise set the scale for everything else.
    xlims!(axis, low, high)
    ylims!(axis, low, high)
    return figure
end

"The temperature ratio of one system: every measurement, and the systematic trend through them."
function figure_temperature_ratio(result)
    figure = Figure(; size = (SINGLE_COLUMN, 0.72 * SINGLE_COLUMN))
    axis = Axis(
        figure[1, 1]; xlabel = L"Heavy fragment mass number $A_H$", ylabel = L"R_T = T_L / T_H"
    )
    hlines!(axis, [1.0]; color = :black, linestyle = :dashdot, linewidth = 0.7)
    for curve in result.R_T
        isempty(curve) && continue
        scatter!(axis, curve.A_H, curve.value; color = (:grey40, 0.45), markersize = 3.2)
    end
    trend = systematic_trend(result).R_T
    band!(
        axis,
        trend.A_H,
        max.(trend.value .- trend.σ, 0.0),
        trend.value .+ trend.σ;
        color = (RGBf(0.835, 0.369, 0.0), 0.25),
    )
    lines!(axis, trend.A_H, trend.value; color = RGBf(0.835, 0.369, 0.0), linewidth = 1.6)
    text!(
        axis,
        0.97,
        0.93;
        text = "$(result.configuration.system.label)\n$(length(result.data_sets)) measurements\nsystematic trend",
        space = :relative,
        align = (:right, :top),
        fontsize = 7,
    )
    ylims!(axis, 0, 3)
    return figure
end

"From the sawtooth to the temperature ratio, one measurement, on a shared abscissa."
function figure_method_chain(result, label)
    index = findfirst(set -> set.label == label, result.data_sets)
    index === nothing && (index = argmax(length.(result.r_ν)))
    data = result.data_sets[index]
    A₀ = result.configuration.system.A₀
    colour = RGBf(0.0, 0.447, 0.698)

    figure = Figure(; size = (SINGLE_COLUMN, 1.25 * SINGLE_COLUMN))
    axes = [
        Axis(figure[row, 1]; ylabel = ylabel, xticklabelsvisible = row == 3) for
        (row, ylabel) in enumerate((L"\nu", L"r_\nu = \nu_H/(\nu_L + \nu_H)", L"R_T = T_L/T_H"))
    ]
    axes[3].xlabel = L"Heavy fragment mass number $A_H$"
    linkxaxes!(axes...)
    # Enough clearance that the lowest tick label of one panel cannot meet the highest of the next.
    rowgap!(figure.layout, 10)

    masses = A_H_range(result.configuration)
    heavy = [(A, ν) for (A, ν) in zip(data.A, data.ν) if A in masses]
    light = [(A₀ - A, ν) for (A, ν) in zip(data.A, data.ν) if A₀ - A in masses && 2A < A₀]
    scatter!(
        axes[1], first.(heavy), last.(heavy); color = colour, markersize = 3.5, label = L"\nu_H"
    )
    scatter!(
        axes[1],
        first.(light),
        last.(light);
        color = RGBf(0.835, 0.369, 0.0),
        markersize = 3.5,
        marker = :rect,
        label = L"\nu_L",
    )
    axislegend(
        axes[1];
        position = :lt,
        framevisible = false,
        labelsize = 6.5,
        orientation = :horizontal,
        patchsize = (8, 6),
    )

    hlines!(axes[2], [0.5]; color = :black, linestyle = :dashdot, linewidth = 0.7)
    scatter!(
        axes[2],
        result.r_ν[index].A_H,
        result.r_ν[index].value;
        color = colour,
        markersize = 3.5,
    )
    hlines!(axes[3], [1.0]; color = :black, linestyle = :dashdot, linewidth = 0.7)
    scatter!(
        axes[3],
        result.R_T[index].A_H,
        result.R_T[index].value;
        color = colour,
        markersize = 3.5,
    )
    text!(
        axes[1],
        0.97,
        0.92;
        text = "$(result.configuration.system.label) · $(data.label)",
        space = :relative,
        align = (:right, :top),
        fontsize = 6.5,
        color = :grey40,
    )
    ylims!(axes[2], 0, 1)
    ylims!(axes[3], 0, 3)
    return figure
end

"How much the level density prescription moves the answer."
function figure_prescriptions(bsfg, gc)
    figure = Figure(; size = (SINGLE_COLUMN, 0.72 * SINGLE_COLUMN))
    axis = Axis(
        figure[1, 1]; xlabel = L"Heavy fragment mass number $A_H$", ylabel = L"R_T = T_L / T_H"
    )
    hlines!(axis, [1.0]; color = :black, linestyle = :dashdot, linewidth = 0.7)
    for (result, name, colour, style) in (
        (bsfg, "Back-shifted Fermi gas", RGBf(0.0, 0.447, 0.698), :solid),
        (gc, "Gilbert-Cameron", RGBf(0.835, 0.369, 0.0), :dash),
    )
        trend = systematic_trend(result).R_T
        lines!(
            axis,
            trend.A_H,
            trend.value;
            color = colour,
            linewidth = 1.6,
            linestyle = style,
            label = name,
        )
    end
    axislegend(axis; position = :rt, framevisible = false, labelsize = 6.5, patchsize = (10, 6))
    text!(
        axis,
        0.03,
        0.05;
        text = "$(bsfg.configuration.system.label)\nsystematic trend under each prescription",
        space = :relative,
        align = (:left, :bottom),
        fontsize = 6.5,
        color = :grey40,
    )
    return figure
end

function main()
    mkpath(ASSETS)
    # The active project here is `docs/`, so DrWatson's default data directory would resolve
    # under it; the data lives beside the package.
    configuration = load_configuration(CONFIG; data_directory = DATA_DIRECTORY)
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

    # The static figures need the pipeline run rather than one fitted curve. Output is not written
    # again here; these read the results in memory.
    results = Dict{String,PipelineResult}()
    for case in ("U233_nf", "U235_nf", "Cf252_0f")
        settings = load_configuration(
            joinpath(dirname(@__DIR__), "config", "$(case).toml");
            data_directory = DATA_DIRECTORY,
        )
        results[case] = run_pipeline(settings; write_output = false)
    end
    # The same system under the other prescription, for the comparison of the two.
    source = joinpath(dirname(@__DIR__), "config", "U233_nf.toml")
    gilbert_cameron = mktempdir() do directory
        path = joinpath(directory, "gilbert_cameron.toml")
        write(
            path,
            replace(read(source, String), "prescription = \"BSFG\"" => "prescription = \"GC\""),
        )
        return run_pipeline(
            load_configuration(path; data_directory = DATA_DIRECTORY); write_output = false
        )
    end

    with_theme(publication_theme()) do
        for (name, figure) in (
            ("published_comparison.png", figure_published_comparison(results)),
            ("temperature_ratio.png", figure_temperature_ratio(results["Cf252_0f"])),
            ("method_chain.png", figure_method_chain(results["Cf252_0f"], "A. Goeoek 2014")),
            (
                "level_density_prescriptions.png",
                figure_prescriptions(results["U233_nf"], gilbert_cameron),
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
