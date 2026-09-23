# The handoff to a consuming code, exercised rather than assumed.
#
# Everything below reads only what a run writes — the manifest and the files it names — using no
# internal function of the package. That is the point: this is the contract a prompt emission code
# depends on, and a change that breaks it would otherwise be found by whoever is downstream rather
# than here. The rules are those of the consumer's manifest reader, restated on the producer's
# side so that the contract is asserted where the manifest is written.

"A consumer's view of one segmented curve: what to call it, and R_T at every mass number."
struct ConsumedCurve
    label::String
    kind::String
    A_H::Vector{Int}
    R_T::Vector{Float64}
    uncertainty::Vector{Float64}
end

# The one manifest a staged run directory holds. A consumer selects it by the `manifest_` prefix
# and refuses a directory holding none or several, so a run directory must carry exactly one.
function staged_manifest(directory::AbstractString)
    names = filter(f -> startswith(f, "manifest_") && endswith(f, ".toml"), readdir(directory))
    length(names) == 1 || error("expected one manifest in $(directory), found $(names)")
    return joinpath(directory, only(names))
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
    path = staged_manifest(directory)
    manifest = TOML.parsefile(path)

    curves = ConsumedCurve[]
    headers = Vector{String}[]
    for entry in manifest["segmented_curve"]
        haskey(entry, "temperature_ratio_file") || continue
        table = CSV.read(joinpath(dirname(path), entry["temperature_ratio_file"]), DataFrame)
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
        system = manifest["system"],
        run = manifest["run"],
        entries = manifest["segmented_curve"],
        curves = curves,
        headers = headers,
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
        result = run_pipeline(configuration)
        # Staged as the consumer stages it: the whole run directory, copied.
        directory = joinpath(root, "sims", "U233_nth", run_identifier(configuration))
        write_results(result, directory)
        run = consume_run(directory)

        @testset "the manifest satisfies every rule of the consumer's reader" begin
            # [system]: every key present with its type; recorded, not matched.
            for (key, type) in (
                ("label", String),
                ("notation", String),
                ("target_A", Int),
                ("target_Z", Int),
                ("channel", String),
                ("reaction", String),
                ("incident_energy_MeV", Real),
                ("compound_A", Int),
                ("compound_Z", Int),
            )
                @test haskey(run.system, key)
                @test run.system[key] isa type
            end
            @test isfinite(run.system["incident_energy_MeV"])
            # [run]: the quantity declared, so r_ν cannot be staged as R_T.
            @test run.run["ordinate"] == "R_T"
            @test run.run["abscissa"] == ["A_H"]
            # [[segmented_curve]]: at least one; label non-empty and unique, since a label is
            # what selects a curve; kind from the closed set; the two files distinct.
            @test !isempty(run.entries)
            labels = [entry["label"] for entry in run.entries]
            @test all(!isempty, labels)
            @test allunique(labels)
            for entry in run.entries
                @test entry["kind"] in ("dataset", "systematic_trend")
                @test entry["temperature_ratio_file"] isa String
                @test entry["multiplicity_ratio_pivots_file"] isa String
                @test entry["temperature_ratio_file"] != entry["multiplicity_ratio_pivots_file"]
                # Paths resolve against the manifest's own directory.
                @test isfile(joinpath(directory, entry["temperature_ratio_file"]))
                @test isfile(joinpath(directory, entry["multiplicity_ratio_pivots_file"]))
                @test !isabspath(entry["temperature_ratio_file"])
            end
            # Additive fields the consumer ignores, stating what each curve rests on.
            for entry in run.entries
                @test entry["first_A_H"] ≤ entry["last_A_H"]
                @test entry["pairs"] ≥ 1
                @test 0 < entry["coverage"] ≤ 1
            end
            # Exactly one manifest per run directory, selected by its prefix.
            @test count(startswith("manifest_"), readdir(directory)) == 1
            @test basename(staged_manifest(directory)) ==
                "manifest_$(run_identifier(configuration)).toml"
        end

        @testset "the system is identified without parsing a label" begin
            @test run.system["label"] == "U233_nth"
            @test run.system["notation"] == "²³³U(nth,f)"
            @test run.system["target_A"] == 233
            @test run.system["target_Z"] == 92
            @test run.system["channel"] == "nth"
            @test run.system["reaction"] == "n,f"
            @test run.system["compound_A"] == 234
            @test run.system["compound_Z"] == 92
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
                # Every uncertainty of a fitted curve is a number: zero at a pin, positive
                # elsewhere. An unquoted one would be an empty field, and there is none here.
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
            entry = only(filter(e -> e["kind"] == "systematic_trend", run.entries))
            @test entry["first_A_H"] == first(produced.A_H)
            @test entry["last_A_H"] == last(produced.A_H)
            @test entry["pairs"] == length(result.consensus_r_ν)
        end

        @testset "excluded and unfitted datasets are visible, not missing" begin
            # A dataset that supports no fit has no segmented curve, and a consumer would
            # otherwise see only an absence. The diagnostics say which datasets were read and
            # what became of them.
            table = CSV.read(joinpath(directory, "dataset_diagnostics.csv"), DataFrame)
            @test nrow(table) == length(result.datasets)
            @test all(in(("true", "false")), string.(table.pooled))
            @test Set(string.(table.dataset)) ⊇
                Set(c.label for c in run.curves if c.kind == "dataset")
            @test all(!isempty, string.(table.segmented_curve))
        end

        @testset "the result files say which quantity they hold" begin
            # A name states the quantity and the abscissa, so a directory listing is readable and
            # an output can be fed back in without translation. The run token is on the directory
            # and on the manifest, not on every file.
            files = readdir(directory)
            @test any(startswith("R_T_vs_A_H_segmented_"), files)
            @test any(startswith("r_nu_vs_A_H_segmented_"), files)
            @test any(startswith("r_nu_vs_A_H_pivots_"), files)
            @test "r_nu_vs_A_H_consensus_systematic_trend.csv" in files
            @test "total_average_R_T.csv" in files
            @test "dataset_diagnostics.csv" in files
            @test !any(contains("parameterized"), files)
            @test count(contains(run_identifier(configuration)), files) == 1
            averages = CSV.read(joinpath(directory, "total_average_R_T.csv"), DataFrame)
            @test names(averages) == [
                "segmented_curve",
                "mass_yield",
                "R_T",
                "R_T_uncertainty",
                "R_T_uncertainty_independent_points",
            ]
            @test all(averages.R_T_uncertainty .≥ 0)
        end
    end
end
