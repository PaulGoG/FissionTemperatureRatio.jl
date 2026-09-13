# The handoff to a consuming code, exercised rather than assumed.
#
# Everything below reads only what a run writes — the manifest and the files it names — using no
# internal function of the package. That is the point: this is the contract a prompt emission code
# depends on, and a change that breaks it would otherwise be found by whoever is downstream rather
# than here.

"A consumer's view of one segmented curve: what to call it, and R_T at every mass number."
struct ConsumedCurve
    label::String
    kind::String
    A_H::Vector{Int}
    R_T::Vector{Float64}
    uncertainty::Vector{Float64}
end

"""
Read a run the way a downstream code would: find the manifest, take the system it describes, and
load the temperature ratio each segmented curve offers.

Deliberately naive, and deliberately positional. It knows the manifest key names and that a ratio
table holds the abscissa, the quantity and its uncertainty in that order — **not** the column
names. Nothing downstream may be coupled to header text, which is what makes a header rename a
no-op for every consumer.
"""
function consume_run(directory::AbstractString)
    manifests = filter(
        f -> startswith(f, "manifest_") && endswith(f, ".toml"), readdir(directory)
    )
    length(manifests) == 1 || error("expected one manifest in $(directory), found $(manifests)")
    manifest = TOML.parsefile(joinpath(directory, first(manifests)))

    curves = ConsumedCurve[]
    headers = Vector{String}[]
    for entry in manifest["segmented_curve"]
        haskey(entry, "temperature_ratio_file") || continue
        table = CSV.read(joinpath(directory, entry["temperature_ratio_file"]), DataFrame)
        push!(headers, String.(names(table)))
        push!(
            curves,
            ConsumedCurve(
                entry["label"],
                entry["kind"],
                Vector{Int}(table[!, 1]),
                Vector{Float64}(table[!, 2]),
                Vector{Float64}(table[!, 3]),
            ),
        )
    end
    return (
        system = manifest["system"], run = manifest["run"], curves = curves, headers = headers
    )
end

"""
Evaluate a consumed curve at a mass number, by the convention a consuming code uses: linear
between tabulated points, undefined outside them.
"""
function evaluate_consumed(curve::ConsumedCurve, A::Integer)
    (A < first(curve.A_H) || A > last(curve.A_H)) && return nothing
    index = findfirst(≥(A), curve.A_H)
    curve.A_H[index] == A && return curve.R_T[index]
    lower = index - 1
    slope = (curve.R_T[index] - curve.R_T[lower]) / (curve.A_H[index] - curve.A_H[lower])
    return curve.R_T[lower] + slope * (A - curve.A_H[lower])
end

DATA_AVAILABLE && @testset "the handoff to a consuming code" begin
    configuration = load_configuration(
        joinpath(pkgdir(FissionTemperatureRatio), "config", "U233_nth.toml");
        data_directory = DATA_DIRECTORY,
    )

    mktempdir() do root
        result = run_pipeline(configuration; output_root = root)
        directory = joinpath(root, "results", configuration.output.subdirectory)
        run = consume_run(directory)

        @testset "the system is identified without parsing a label" begin
            @test run.system["label"] == "U233_nth"
            @test run.system["notation"] == "²³³U(nth,f)"
            @test run.system["target_A"] == 233
            @test run.system["target_Z"] == 92
            @test run.system["channel"] == "nth"
            @test run.system["reaction"] == "n,f"
            @test run.system["compound_A"] == 234
            @test run.system["compound_Z"] == 92
            @test run.run["ordinate"] == "R_T"
            @test run.run["abscissa"] == ["A_H"]
        end

        @testset "a column is named for its quantity, and read by position" begin
            # The manifest states the columns so that a reader knows what the file holds; the
            # reader above took them by position, and the two must agree.
            @test run.run["columns"] == ["A_H", "R_T", "R_T_uncertainty"]
            @test !isempty(run.headers)
            for header in run.headers
                @test header == run.run["columns"]
                @test "value" ∉ header
                @test "uncertainty" ∉ header
            end
        end

        @testset "every segmented curve names a file that exists and parses" begin
            @test !isempty(run.curves)
            @test count(c -> c.kind == "systematic_trend", run.curves) == 1
            for curve in run.curves
                @test !isempty(curve.A_H)
                @test length(curve.A_H) == length(curve.R_T) == length(curve.uncertainty)
                @test issorted(curve.A_H)
                @test allunique(curve.A_H)
                @test all(>(0), curve.R_T)
                @test all(≥(0), curve.uncertainty)
            end
        end

        @testset "the tabulation is dense, so interpolation is exact" begin
            # The reason the manifest points a consumer at this file rather than at the segment
            # pivots: R_T is not piecewise-linear even where the multiplicity ratio is, so
            # interpolating between the pivots would cut across the structure the level density
            # parameter ratio puts into it. Tabulated at every mass number, any interpolation a
            # consumer uses reproduces the tabulated value exactly.
            for curve in run.curves
                @test curve.A_H == collect(first(curve.A_H):last(curve.A_H))
                for (index, A) in enumerate(curve.A_H)
                    @test evaluate_consumed(curve, A) == curve.R_T[index]
                end
                @test evaluate_consumed(curve, last(curve.A_H) + 1) === nothing
                @test evaluate_consumed(curve, first(curve.A_H) - 1) === nothing
            end
        end

        @testset "the exact identity at the symmetric split survives the round trip" begin
            # R_T(A₀/2) = 1 holds by construction; a consumer must see it in the file, not be
            # told about it.
            for curve in run.curves
                index = findfirst(==(configuration.system.A₀ ÷ 2), curve.A_H)
                index === nothing && continue
                @test curve.R_T[index] ≈ 1 atol = 1e-10
            end
        end

        @testset "what the manifest says matches what the run produced" begin
            trend = only(filter(c -> c.kind == "systematic_trend", run.curves))
            @test trend.label == TREND_LABEL
            produced = systematic_trend(result).R_T
            @test trend.A_H == produced.A_H
            # The file is rounded to the configured number of significant figures; the comparison
            # honours that rather than demanding the file carry more precision than it claims.
            # Rounding to n significant figures moves a value by at most 5·10⁻ⁿ of itself.
            tolerance = 5 * 10.0^(-configuration.output.significant_digits)
            @test all(isapprox.(trend.R_T, produced.ratio; rtol = tolerance))
        end

        @testset "excluded and unfitted datasets are visible, not missing" begin
            # A dataset that supports no fit has no segmented curve, and a consumer would
            # otherwise see only an absence. The diagnostics say which datasets were read and
            # what became of them.
            diagnostics = only(
                filter(f -> startswith(f, "dataset_diagnostics_"), readdir(directory))
            )
            table = CSV.read(joinpath(directory, diagnostics), DataFrame)
            @test nrow(table) == length(result.datasets)
            @test all(in(("true", "false")), string.(table.pooled))
            @test Set(string.(table.dataset)) ⊇
                Set(c.label for c in run.curves if c.kind == "dataset")
        end

        @testset "the result files say which quantity they hold" begin
            # A name states the quantity and the abscissa, so a directory listing is readable and
            # an output can be fed back in without translation.
            files = readdir(directory)
            @test any(startswith("R_T_vs_A_H_segmented_"), files)
            @test any(startswith("r_nu_vs_A_H_segmented_"), files)
            @test any(startswith("r_nu_vs_A_H_pivots_"), files)
            @test "r_nu_vs_A_H_consensus_systematic_trend.csv" in files
            @test any(startswith("total_average_R_T_"), files)
            @test !any(contains("parameterized"), files)
        end
    end
end
