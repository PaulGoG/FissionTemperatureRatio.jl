# How the datasets of one experiment are combined before they are pooled: alternative analyses of
# the same events by their mean, a republication through its successor alone, anything else as a
# pool is combined. Synthetic input only, so this runs on a bare clone.

const EXPERIMENT_FILES = (
    "20000001_A.First_1979.dat", "20000002_B.First_1979.dat", "20000003_C.Other_2000.dat"
)

# ν(A) of both fragments for the heavy masses 127:140 and the symmetric split, from the reference
# ratio displaced by a constant and a total of four neutrons per fission.
function write_displaced_multiplicity(path::AbstractString, offset::Real)
    rows = Dict{Int,Float64}(126 => 2.0)
    for A_H in 127:140
        r = reference_ratio(A_H) + offset + 0.002 * iseven(A_H)
        rows[A_H] = 4 * r
        rows[252 - A_H] = 4 * (1 - r)
    end
    open(path, "w") do io
        println(io, "A nu nu_uncertainty")
        for A in sort!(collect(keys(rows)))
            println(io, A, " ", rows[A], " 0.05")
        end
    end
    return path
end

# The record lines naming the partners of a dataset and, where given, its relation to them.
function correlation_lines(partners::AbstractString...; relation::AbstractString = "")
    quoted = join(("\"$(p)\"" for p in partners), ", ")
    lines = ["correlated_with = [$(quoted)]"]
    isempty(relation) || push!(lines, "correlation_relation = \"$(relation)\"")
    return join(lines, "\n")
end

# The record lines of a dataset whose 20 non-integer masses were interpolated onto `rows` rows.
function interpolated_lines(rows::Integer)
    return "mass_treatment = \"interpolated\"\nmass_values_non_integer = 20\n" *
           "rows_written = $(rows)"
end

const ALTERNATIVE_ANALYSES = (
    correlation_lines("20000002"; relation = "alternative_analysis"),
    correlation_lines("20000001"; relation = "alternative_analysis"),
    "",
)

const ONE_EXPERIMENT = [["A. First 1979", "B. First 1979"]]

# A run over the datasets of `EXPERIMENT_FILES`, the reference ratio displaced by +0.02, +0.03 and
# -0.03, whose retrieval record carries, for each file in turn, the lines of `entries` and the
# qualifiers of `qualifiers`. `excluded` is an accession the configuration keeps out of the pooling.
function experiment_run(
    directory::AbstractString,
    entries::NTuple{3,String};
    qualifiers::NTuple{3,String} = ("[]", "[]", "[]"),
    excluded::Union{Nothing,String} = nothing,
    leave_one_out::Bool = false,
)
    # Vanishing shell corrections, as in the pipeline tests. The reader drops trailing zeros of
    # S(Z) as padding, so the last row, at a nucleon number no fragment carries, is not zero.
    open(joinpath(directory, "shell_corrections.dat"), "w") do io
        println(io, "n S_N S_Z")
        for n in 20:119
            println(io, n, " 0.0 0.0")
        end
        return println(io, 120, " 1.0 1.0")
    end
    held = joinpath(directory, "experiment")
    mkpath(held)
    record = IOBuffer()
    offsets = (0.02, 0.03, -0.03)
    for (file, offset, entry, tags) in zip(EXPERIMENT_FILES, offsets, entries, qualifiers)
        write_displaced_multiplicity(joinpath(held, file), offset)
        println(record, "[[accepted]]")
        println(record, "file = \"$(file)\"")
        println(record, "identifier = \"$(first(split(file, '_')))\"")
        println(record, "qualifiers = $(tags)")
        isempty(entry) || println(record, entry)
        println(record)
    end
    println(record, "[run]")
    println(record, "package_version = \"0.2.7\"")
    write(joinpath(held, "retrieval.toml"), take!(record))
    multiplicity = "subdirectory = \"experiment\""
    if excluded !== nothing
        multiplicity *= "\nexclude = [{ accession = \"$(excluded)\", reason = \"test\" }]"
    end
    path = joinpath(directory, "configuration.toml")
    write(path, replace(PIPELINE_CONFIGURATION, "subdirectory = \"datasets\"" => multiplicity))
    return run_pipeline(
        load_configuration(path; data_directory = directory); leave_one_out = leave_one_out
    )
end

# Grouped, and combined by the rule of a pool: the members first, then their combination with
# the third dataset. With pooling weights of one the standard error of the members' combination
# is its uncertainty as one measurement, so the composition through `consensus` is exact.
function test_combined_as_a_pool(result)
    first_, second, other = result.r_ν
    @test result.correlation_groups == ONE_EXPERIMENT
    @test result.consensus_r_ν.ratio ≈ consensus([consensus([first_, second]), other]).ratio
    return nothing
end

function test_alternative_analyses_pooled()
    @testset "alternative analyses stated by every member enter the pool as their mean" begin
        mktempdir() do directory
            result = experiment_run(directory, ALTERNATIVE_ANALYSES)
            first_, second, other = result.r_ν
            @test result.correlation_groups == ONE_EXPERIMENT
            either = first(
                FissionTemperatureRatio._alternative_analyses(
                    [first_, second], [1.0, 1.0], "either"
                ),
            )
            @test result.consensus_r_ν.ratio ≈ consensus([either, other]).ratio
            # Pooled side by side, the two analyses would weigh as two measurements.
            @test !isapprox(
                result.consensus_r_ν.ratio[2:end],
                consensus([first_, second, other]).ratio[2:end];
                rtol = 1e-6,
            )
        end
    end
    return nothing
end

function test_relations_not_shared()
    @testset "a relation stated by one member alone is not acted on" begin
        mktempdir() do directory
            result = experiment_run(
                directory,
                (
                    correlation_lines("20000002"; relation = "alternative_analysis"),
                    correlation_lines("20000001"),
                    "",
                ),
            )
            test_combined_as_a_pool(result)
        end
    end
    @testset "members stating different relations are combined as a pool" begin
        mktempdir() do directory
            result = experiment_run(
                directory,
                (
                    correlation_lines("20000002"; relation = "alternative_analysis"),
                    correlation_lines("20000001"; relation = "repeated_run"),
                    "",
                ),
            )
            test_combined_as_a_pool(result)
        end
    end
    @testset "members stating no relation are combined as a pool" begin
        mktempdir() do directory
            result = experiment_run(
                directory, (correlation_lines("20000002"), correlation_lines("20000001"), "")
            )
            test_combined_as_a_pool(result)
        end
    end
    return nothing
end

function test_partners_outside_the_pool()
    @testset "a partner the directory does not hold forms no experiment" begin
        mktempdir() do directory
            result = experiment_run(directory, (correlation_lines("29999999"), "", ""))
            first_, second, other = result.r_ν
            @test isempty(result.correlation_groups)
            @test result.consensus_r_ν.ratio ≈ consensus([first_, second, other]).ratio
        end
    end
    @testset "a partner excluded by the configuration forms no experiment" begin
        mktempdir() do directory
            result = experiment_run(directory, ALTERNATIVE_ANALYSES; excluded = "20000002")
            first_, _, other = result.r_ν
            @test isempty(result.correlation_groups)
            @test pooled_datasets(result) == ["A. First 1979", "C. Other 2000"]
            @test result.consensus_r_ν.ratio ≈ consensus([first_, other]).ratio
        end
    end
    return nothing
end

function test_republication()
    republished = (
        correlation_lines("20000002"; relation = "republication"),
        correlation_lines("20000001"; relation = "republication"),
        "",
    )
    marked = "[\"superseded: by 20000002\"]"
    @testset "a republication enters the pool through its one unmarked member" begin
        mktempdir() do directory
            result = experiment_run(directory, republished; qualifiers = (marked, "[]", "[]"))
            _, second, other = result.r_ν
            @test superseded_datasets(result) == Dict("A. First 1979" => "B. First 1979")
            @test pooled_datasets(result) == ["B. First 1979", "C. Other 2000"]
            @test result.consensus_r_ν.ratio ≈ consensus([second, other]).ratio
        end
    end
    @testset "a republication with two unmarked members supersedes nothing" begin
        mktempdir() do directory
            result = experiment_run(directory, republished)
            @test isempty(superseded_datasets(result))
            test_combined_as_a_pool(result)
        end
    end
    @testset "a republication with no unmarked member supersedes nothing" begin
        mktempdir() do directory
            result = experiment_run(directory, republished; qualifiers = (marked, marked, "[]"))
            @test isempty(superseded_datasets(result))
            test_combined_as_a_pool(result)
        end
    end
    @testset "a dataset superseded by an excluded successor is pooled with neither" begin
        mktempdir() do directory
            warning = (:warn, r"superseded by a republication that is not pooled")
            result = @test_logs warning match_mode = :any experiment_run(
                directory, republished; qualifiers = (marked, "[]", "[]"), excluded = "20000002"
            )
            @test superseded_datasets(result) == Dict("A. First 1979" => "B. First 1979")
            @test pooled_datasets(result) == ["C. Other 2000"]
        end
    end
    return nothing
end

function test_interpolated_alternative_analyses()
    @testset "interpolated analyses enter at the mean of their pooling factors" begin
        mktempdir() do directory
            # Pooling weights 20/29 and 1/2.
            entries = (
                ALTERNATIVE_ANALYSES[1] * "\n" * interpolated_lines(29),
                ALTERNATIVE_ANALYSES[2] * "\n" * interpolated_lines(40),
                "",
            )
            result = experiment_run(directory, entries; leave_one_out = true)
            first_, second, other = result.r_ν
            @test result.correlation_groups == ONE_EXPERIMENT
            @test [entry.datasets for entry in result.leave_one_out] == [["A. First 1979", "B. First 1979"], ["C. Other 2000"]]
            either, measured = FissionTemperatureRatio._alternative_analyses(
                [first_, second], [20 / 29, 0.5], "either"
            )
            @test all(measured .≈ (20 / 29 + 0.5) / 2)
            # The two-unit pool, the experiment carrying its measured fraction point by point.
            expected = first(
                FissionTemperatureRatio._consensus([either, other], "expected", [measured, 1.0])
            )
            @test result.consensus_r_ν.ratio ≈ expected.ratio
            @test result.consensus_r_ν.σ[2:end] ≈ expected.σ[2:end]
        end
    end
    return nothing
end

@testset "experiments in the pool" begin
    a = RatioCurve([130, 131], [0.30, 0.40], [0.010, 0.020], "a")
    b = RatioCurve([130, 131], [0.34, 0.40], [0.030, 0.010], "b")
    a₃ = RatioCurve([130, 131, 132], [0.30, 0.40, 0.50], [0.01, 0.02, 0.03], "a")
    analyses(curves, factors) =
        FissionTemperatureRatio._alternative_analyses(curves, factors, "either")

    @testset "two analyses are averaged, at the larger uncertainty and half their difference" begin
        curve, measured = analyses([a, b], [1.0, 1.0])
        @test curve.A_H == [130, 131]
        @test curve.ratio ≈ [0.32, 0.40]
        @test curve.σ[1] ≈ hypot(0.030, 0.02)
        @test curve.σ[2] ≈ 0.020
        @test measured == [1.0, 1.0]
        @test curve.label == "either"
    end

    @testset "a mass number one analysis alone holds keeps its value and uncertainty" begin
        curve, measured = analyses([a₃, b], [1.0, 1.0])
        @test curve.A_H == [130, 131, 132]
        @test curve.ratio[3] == 0.50
        @test curve.σ[3] == 0.03
        @test measured[3] == 1.0
    end

    @testset "a point quotes an uncertainty where an analysis does, none where none does" begin
        quoting = RatioCurve([130, 131], [0.30, 0.30], [missing, 0.02], "quoting")
        silent = RatioCurve(
            [130, 131], [0.36, 0.36], Union{Missing,Float64}[missing, missing], "silent"
        )
        curve, _ = analyses([quoting, silent], [1.0, 1.0])
        # The values differ, and still no uncertainty is invented from their difference alone.
        @test ismissing(curve.σ[1])
        @test curve.ratio[1] ≈ 0.33
        @test curve.σ[2] ≈ hypot(0.02, 0.03)
    end

    @testset "of three analyses the half difference is that of the extremes" begin
        curves = [
            RatioCurve([130], [value], [σ], "member") for
            (value, σ) in ((0.30, 0.01), (0.33, 0.04), (0.39, 0.02))
        ]
        curve, _ = analyses(curves, [1.0, 1.0, 1.0])
        @test curve.ratio[1] ≈ 0.34
        @test curve.σ[1] ≈ hypot(0.04, 0.045)
    end

    @testset "the measured fraction is the mean pooling factor of the analyses at a point" begin
        _, measured = analyses([a₃, b], [1.0, 0.5])
        @test measured ≈ [0.75, 0.75, 1.0]
        _, measured = analyses([a₃, b], [0.5, 1.0])
        @test measured ≈ [0.75, 0.75, 0.5]
    end

    @testset "the uncertainty is not reduced as that of independent measurements is" begin
        curve, _ = analyses([a, b], [1.0, 1.0])
        pool = consensus([a, b])
        @test curve.σ[1] > pool.σ[1]
        # Equal values: the larger quoted uncertainty, where a pool would go below the smaller.
        @test curve.σ[2] ≈ 0.020
        @test pool.σ[2] < 0.010
        @test curve.σ[2] > pool.σ[2]
    end

    test_alternative_analyses_pooled()
    test_relations_not_shared()
    test_partners_outside_the_pool()
    test_republication()
    test_interpolated_alternative_analyses()
end
