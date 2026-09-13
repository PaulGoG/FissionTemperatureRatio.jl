# The handoff to a consuming code, exercised rather than assumed.
#
# Everything below reads only what a run writes — the manifest and the files it names — using no
# internal function of the package. That is the point: this is the contract a prompt emission code
# depends on, and a change that breaks it would otherwise be found by whoever is downstream rather
# than here.

"A consumer's view of one parameterization: what to call it, and R_T at every mass number."
struct ConsumedCurve
    label::String
    kind::String
    A_H::Vector{Int}
    R_T::Vector{Float64}
    uncertainty::Vector{Float64}
end

"""
Read a run the way a downstream code would: find the manifest, take the system it describes, and
load the temperature ratio each parameterization offers.

Deliberately naive. It knows the manifest key names and the three column names, nothing else.
"""
function consume_run(directory::AbstractString)
    manifests = filter(
        f -> startswith(f, "manifest_") && endswith(f, ".toml"), readdir(directory)
    )
    length(manifests) == 1 || error("expected one manifest in $(directory), found $(manifests)")
    manifest = TOML.parsefile(joinpath(directory, first(manifests)))

    curves = ConsumedCurve[]
    for entry in manifest["parameterization"]
        haskey(entry, "temperature_ratio_file") || continue
        table = CSV.read(joinpath(directory, entry["temperature_ratio_file"]), DataFrame)
        push!(
            curves,
            ConsumedCurve(
                entry["label"],
                entry["kind"],
                Vector{Int}(table.A_H),
                Vector{Float64}(table.value),
                Vector{Float64}(table.uncertainty),
            ),
        )
    end
    return (system = manifest["system"], run = manifest["run"], curves = curves)
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
        joinpath(pkgdir(FissionTemperatureRatio), "config", "U233_nf.toml");
        data_directory = DATA_DIRECTORY,
    )

    mktempdir() do root
        result = run_pipeline(configuration; output_root = root)
        directory = joinpath(root, "results", configuration.output.subdirectory)
        run = consume_run(directory)

        @testset "the system is identified without parsing a label" begin
            @test run.system["label"] == "U233_nf"
            @test run.system["target_A"] == 233
            @test run.system["target_Z"] == 92
            @test run.system["reaction"] == "n,f"
            @test run.system["compound_A"] == 234
            @test run.system["compound_Z"] == 92
            @test run.run["columns"] == ["A_H", "value", "uncertainty"]
        end

        @testset "every parameterization names a file that exists and parses" begin
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
            # The file is rounded to the configured number of digits; the comparison honours that
            # rather than demanding the file carry more precision than it claims.
            tolerance = 10.0^(-configuration.output.digits) / 2
            @test all(abs.(trend.R_T .- produced.value) .≤ tolerance)
        end

        @testset "excluded and unfitted sets are visible, not missing" begin
            # A set that supports no fit has no parameterization, and a consumer would otherwise
            # see only an absence. The diagnostics say which sets were read and what became of
            # them.
            diagnostics = only(filter(f -> startswith(f, "diagnostics_"), readdir(directory)))
            table = CSV.read(joinpath(directory, diagnostics), DataFrame)
            @test nrow(table) == length(result.data_sets)
            @test all(in(("true", "false")), string.(table.pooled))
            @test Set(string.(table.set)) ⊇
                Set(c.label for c in run.curves if c.kind == "data_set")
        end
    end
end
