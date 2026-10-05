# The pipeline on synthetic input: no shipped data is needed, so this runs on a bare clone.

const PIPELINE_CONFIGURATION = """
[system]
target_A = 252
target_Z = 98
channel = "sf"

[fragmentation]
charges_per_mass = 5
heavy_mass_max = 140

[level_density]
model = "GC"
shell_correction_file = "shell_corrections.dat"
deformed_branch = false

[multiplicity]
subdirectory = "datasets"

[segments]
max_segments = 3
pin_symmetric_split = true

[output]
significant_digits = 6
"""

# ν(A) of both fragments for the heavy masses given, from the reference ratio and a total of four
# neutrons per fission, with a deterministic perturbation so that no fit is exact.
function write_multiplicity(path::AbstractString, heavy_masses; A₀ = 252)
    rows = Dict{Int,Float64}()
    for A_H in heavy_masses
        r = reference_ratio(A_H) + 0.002 * iseven(A_H)
        rows[A_H] = 4 * r
        rows[A₀ - A_H] = 4 * (1 - r)
    end
    open(path, "w") do io
        println(io, "A nu nu_uncertainty")
        for A in sort!(collect(keys(rows)))
            println(io, A, " ", rows[A], " 0.05")
        end
    end
    return path
end

# The testsets below are functions of their own: as part of the one closure that holds the rest
# of this file they made it too long to compile in reasonable time.
function test_trend_covariance(directory, configuration, result, written)
    @testset "the autocorrelation of the datasets' deviations enters the trend's covariance" begin
        # Two datasets that depart from their common mean in opposite directions: point by
        # point in one pair, by a slow drift in the other.
        function write_departing(path, departure)
            rows = Dict{Int,Float64}()
            for A_H in 127:140
                r = reference_ratio(A_H) + departure(A_H)
                rows[A_H] = 4 * r
                rows[252 - A_H] = 4 * (1 - r)
            end
            rows[126] = 2.0
            open(path, "w") do io
                println(io, "A nu nu_uncertainty")
                for A in sort!(collect(keys(rows)))
                    println(io, A, " ", rows[A], " 0.05")
                end
            end
        end
        function autocorrelation(name, departure)
            directory_ = joinpath(directory, name)
            mkpath(directory_)
            write_departing(joinpath(directory_, "up.dat"), A -> departure(A))
            write_departing(joinpath(directory_, "down.dat"), A -> -departure(A))
            path = joinpath(directory, "$(name).toml")
            write(
                path,
                replace(
                    PIPELINE_CONFIGURATION,
                    "subdirectory = \"datasets\"" => "subdirectory = \"$(name)\"",
                ),
            )
            return run_pipeline(load_configuration(path; data_directory = directory))
        end
        alternating = autocorrelation("alternating", A -> 0.01 * (-1)^A)
        drifting = autocorrelation("drifting", A -> 0.002 * (A - 133))
        fit_of(run) = systematic_trend(run).fit

        # Point by point: a negative first lag, taken as no correlation.
        @test first(deviation_correlogram(alternating)) < -0.8
        @test deviation_autocorrelation(alternating) == 0.0
        @test !fit_of(alternating).correlated
        @test fit_of(alternating).expected_wrss == fit_of(alternating).dof

        # A slow drift: the decay fitted to the first four lags enters the trend's covariance.
        ρ = deviation_autocorrelation(drifting)
        @test first(deviation_correlogram(drifting)) > 0.6
        @test ρ == autocorrelation_decay(deviation_correlogram(drifting); lags = 4)
        @test ρ > 0.6
        @test drifting.autocorrelation == ρ
        @test fit_of(drifting).correlated
        @test fit_of(drifting).expected_wrss < fit_of(drifting).dof
        # The values are those of the fit that takes the points as independent; the
        # uncertainty is larger at every mass number but the pinned one.
        independent = FissionTemperatureRatio._fit_trend(
            drifting.r_ν,
            ones(length(drifting.r_ν)),
            (drifting.averaging, drifting.model, drifting.domain),
            drifting.configuration.segments,
            252;
            estimate_autocorrelation = false,
        ).trend
        trend = systematic_trend(drifting)
        @test trend.R_T.ratio == independent.R_T.ratio
        @test trend.fit.breakpoints == independent.fit.breakpoints
        @test trend.R_T.σ[1] ≈ 0 atol = 1e-12
        @test all(trend.R_T.σ[2:end] .> independent.R_T.σ[2:end])
        # The curve of a dataset is fitted to points measured one by one.
        @test all(
            !c.fit.correlated && c.fit.expected_wrss == c.fit.dof for
            c in drifting.segmented_curves if c.kind == "dataset"
        )

        # Identical datasets do not deviate: nothing to estimate, and independent points.
        @test all(ismissing, deviation_correlogram(result))
        @test ismissing(deviation_autocorrelation(result))
        @test !fit_of(result).correlated

        metadata = run_metadata(drifting)["result"]
        record = metadata["trend_uncertainty"]
        @test record["deviation_autocorrelation"] == ρ
        @test record["autocorrelation_estimated"]
        @test record["autocorrelation_lags"] == 4
        @test length(record["deviation_correlogram"]) == CORRELOGRAM_LAGS
        @test record["deviation_correlogram"][1] == first(deviation_correlogram(drifting))
        @test record["reduced_chi_squared"] == trend.fit.wrss / trend.fit.dof
        @test record["expected_wrss"] == trend.fit.expected_wrss
        @test record["chi_squared_over_expectation"] == trend.fit.wrss / trend.fit.expected_wrss
        @test record["covariance_scale"] == max(1, record["chi_squared_over_expectation"])
        @test !haskey(record, "understated_by_about")
        @test metadata["segmented_curves"][SYSTEMATIC_TREND_LABEL]["correlated_points"]
        @test !metadata["segmented_curves"]["up"]["correlated_points"]
        @test run_metadata(drifting)["configuration"]["autocorrelation_lags"] == 4
        # At its default the number of lags is no token of the identifier.
        @test !haskey(run_metadata(drifting)["identifier"]["tokens"], "lags")
        unestimated = run_metadata(result)["result"]["trend_uncertainty"]
        # One type per key: a flag for whether a coefficient could be fitted, and the number
        # the covariance was formed with.
        @test !unestimated["autocorrelation_estimated"]
        @test unestimated["deviation_autocorrelation"] == 0.0
        @test all(isnan, unestimated["deviation_correlogram"])

        table = CSV.read(
            write_results(drifting, joinpath(directory, "output", "drifting"))["segmented_curves"],
            DataFrame,
        )
        row = only(filter(r -> r.label == SYSTEMATIC_TREND_LABEL, eachrow(table)))
        @test row.deviation_autocorrelation ≈ ρ rtol = 1e-5
        @test row.chi_squared_over_expectation ≈ trend.fit.wrss / trend.fit.expected_wrss rtol =
            1e-5
        @test row.chi_squared_over_expectation > row.reduced_chi_squared
        datasets_ = filter(r -> r.kind == "dataset", table)
        @test all(ismissing, datasets_.deviation_autocorrelation)
        @test datasets_.chi_squared_over_expectation == datasets_.reduced_chi_squared

        # The number of lags is a key of the configuration, bounded by the correlogram.
        one_lag = joinpath(directory, "one_lag.toml")
        write(
            one_lag,
            replace(
                PIPELINE_CONFIGURATION,
                "subdirectory = \"datasets\"" => "subdirectory = \"drifting\"",
                "max_segments = 3" => "max_segments = 3\nautocorrelation_lags = 1",
            ),
        )
        first_lag = run_pipeline(load_configuration(one_lag; data_directory = directory))
        @test deviation_autocorrelation(first_lag) ≈ first(deviation_correlogram(drifting)) atol =
            1e-6
        @test run_parameters(first_lag.configuration)["lags"] == 1
        @test run_identifier(first_lag.configuration) != run_identifier(drifting.configuration)
        write(
            one_lag,
            replace(
                PIPELINE_CONFIGURATION,
                "max_segments = 3" => "max_segments = 3\nautocorrelation_lags = 9",
            ),
        )
        @test_throws "segments.autocorrelation_lags" load_configuration(
            one_lag; data_directory = directory
        )
    end
    return nothing
end

function test_leave_one_out(directory, configuration, result, written)
    @testset "the trend is refitted with each pooled dataset left out" begin
        @test [only(entry.datasets) for entry in result.leave_one_out] == pooled_datasets(result)
        @test all(entry.outcome == "segmented curve" for entry in result.leave_one_out)
        @test all(entry.segments ≥ 1 for entry in result.leave_one_out)
        # No yield distribution, so no total average and no spread of one.
        @test all(isempty(entry.total_average_R_T) for entry in result.leave_one_out)
        @test leave_one_out_spread(result, "any") === nothing
        table = CSV.read(written["leave_one_out"], DataFrame)
        @test table.dataset_left_out == pooled_datasets(result)
        @test all(ismissing, table.mass_yield)
        @test all(ismissing, table.R_T)
        @test table.segments == [entry.segments for entry in result.leave_one_out]
        @test basename(written["leave_one_out"]) ==
            "leave_one_out_$(run_identifier(configuration)).csv"
        record = run_metadata(result)["result"]["leave_one_out"]
        @test [only(entry["datasets"]) for entry in record] == pooled_datasets(result)
        @test record[1]["breakpoints"] == result.leave_one_out[1].breakpoints

        # Skipped on request: the curves are the same, and no table is written.
        without = run_pipeline(configuration; leave_one_out = false)
        @test isempty(without.leave_one_out)
        @test systematic_trend(without).R_T.ratio == systematic_trend(result).R_T.ratio
        @test systematic_trend(without).R_T.σ == systematic_trend(result).R_T.σ
        @test !haskey(
            write_results(without, joinpath(directory, "output", "without")), "leave_one_out"
        )
    end
    return nothing
end

function test_experiments_pooled_as_one(directory, configuration, result, written)
    @testset "the datasets of one experiment are pooled as one" begin
        # Two analyses of one experiment and an independent measurement, each departing from
        # the reference by a constant: side by side the experiment would weigh twice.
        experiments = joinpath(directory, "experiments")
        mkpath(experiments)
        function write_offset(name, offset)
            rows = Dict{Int,Float64}(126 => 2.0)
            for A_H in 127:140
                r = reference_ratio(A_H) + offset + 0.002 * iseven(A_H)
                rows[A_H] = 4 * r
                rows[252 - A_H] = 4 * (1 - r)
            end
            open(joinpath(experiments, name), "w") do io
                println(io, "A nu nu_uncertainty")
                for A in sort!(collect(keys(rows)))
                    println(io, A, " ", rows[A], " 0.05")
                end
            end
        end
        write_offset("20000001_A.First_1979.dat", 0.02)
        write_offset("20000002_B.First_1979.dat", 0.03)
        write_offset("20000003_C.Other_2000.dat", -0.03)
        # The relation is named from one side only, and holds for both.
        record(correlated; qualifiers = "[]", second = "") = """
            [[accepted]]
            file = "20000001_A.First_1979.dat"
            identifier = "20000001"
            qualifiers = $(qualifiers)
            $(correlated)

            [[accepted]]
            file = "20000002_B.First_1979.dat"
            identifier = "20000002"
            qualifiers = []
            $(second)

            [[accepted]]
            file = "20000003_C.Other_2000.dat"
            identifier = "20000003"
            qualifiers = []

            [run]
            package_version = "0.2.7"
            """
        path = joinpath(directory, "experiments.toml")
        write(
            path,
            replace(
                PIPELINE_CONFIGURATION,
                "subdirectory = \"datasets\"" => "subdirectory = \"experiments\"",
            ),
        )
        write(joinpath(experiments, "retrieval.toml"), record(""))
        apart = run_pipeline(load_configuration(path; data_directory = directory))
        write(
            joinpath(experiments, "retrieval.toml"),
            record("correlated_with = [\"20000002\"]\ncorrelation_relation = \"repeated_run\""),
        )
        together = run_pipeline(load_configuration(path; data_directory = directory))
        @test correlated_datasets(
            retrieval_record(joinpath(experiments, "20000001_A.First_1979.dat"))
        ) == ["20000002"]
        @test isempty(
            correlated_datasets(
                retrieval_record(joinpath(experiments, "20000003_C.Other_2000.dat"))
            ),
        )

        first_, second, other = together.r_ν
        @test isempty(apart.correlation_groups)
        @test together.correlation_groups == [["A. First 1979", "B. First 1979"]]
        # All three are pooled, and each still offers its own curve.
        @test pooled_datasets(together) == ["A. First 1979", "B. First 1979", "C. Other 2000"]
        @test count(c -> c.kind == "dataset", together.segmented_curves) == 3
        # Side by side: the combination of three. As one: the two of the experiment first,
        # then their combination with the third as two measurements.
        @test apart.consensus_r_ν.ratio ≈ consensus([first_, second, other]).ratio
        experiment = consensus([first_, second])
        expected = consensus([experiment, other])
        @test together.consensus_r_ν.ratio ≈ expected.ratio
        @test together.consensus_r_ν.σ[2:end] ≈ expected.σ[2:end]
        # The independent measurement weighs as much as the experiment, not half as much.
        A = 2:length(other)
        @test all(
            abs.(together.consensus_r_ν.ratio[A] .- other.ratio[A]) .<
            abs.(apart.consensus_r_ν.ratio[A] .- other.ratio[A]),
        )
        @test together.consensus_r_ν.ratio[A] ≈ (experiment.ratio[A] .+ other.ratio[A]) ./ 2 atol =
            2e-3

        # Left out as one, and written as one.
        @test [entry.datasets for entry in together.leave_one_out] == [["A. First 1979", "B. First 1979"], ["C. Other 2000"]]
        @test length(apart.leave_one_out) == 3
        written_ = write_results(together, joinpath(directory, "output", "experiments"))
        refits = CSV.read(written_["leave_one_out"], DataFrame)
        @test refits.dataset_left_out == ["A. First 1979 + B. First 1979", "C. Other 2000"]
        @test string.(refits.accession) == ["20000001 20000002", "20000003"]
        diagnostics = CSV.read(written_["dataset_diagnostics"], DataFrame)
        row(label) = only(filter(r -> r.dataset == label, eachrow(diagnostics)))
        @test row("A. First 1979").pooled_with == "B. First 1979"
        @test row("B. First 1979").pooled_with == "A. First 1979"
        @test ismissing(row("C. Other 2000").pooled_with)
        metadata = run_metadata(together)["result"]
        group = only(metadata["correlation_groups"])
        @test group["datasets"] == ["A. First 1979", "B. First 1979"]
        @test group["accessions"] == ["20000001", "20000002"]
        @test group["relation"] == "repeated_run"
        @test metadata["leave_one_out"][1]["accessions"] == ["20000001", "20000002"]
        @test isempty(run_metadata(apart)["result"]["correlation_groups"])

        # Two reductions of the same events share their statistical errors: the mean, at the
        # larger uncertainty with half their difference added in quadrature.
        write(
            joinpath(experiments, "retrieval.toml"),
            record(
                "correlated_with = [\"20000002\"]\ncorrelation_relation = \"alternative_analysis\"";
                second = "correlated_with = [\"20000001\"]\ncorrelation_relation = \"alternative_analysis\"",
            ),
        )
        reductions = run_pipeline(load_configuration(path; data_directory = directory))
        @test reductions.correlation_groups == [["A. First 1979", "B. First 1979"]]
        @test only(run_metadata(reductions)["result"]["correlation_groups"])["relation"] ==
            "alternative_analysis"
        halves = abs.(first_.ratio .- second.ratio) ./ 2
        either = RatioCurve(
            first_.A_H,
            (first_.ratio .+ second.ratio) ./ 2,
            hypot.(max.(first_.σ, second.σ), halves),
            "either",
        )
        @test reductions.consensus_r_ν.ratio ≈ consensus([either, other]).ratio
        @test reductions.consensus_r_ν.σ[2:end] ≈ consensus([either, other]).σ[2:end]
        # Wider than the inverse-variance combination of two independent values.
        @test all(either.σ[A] .> experiment.σ[A])
        @test [entry.datasets for entry in reductions.leave_one_out] == [["A. First 1979", "B. First 1979"], ["C. Other 2000"]]

        # One result published twice: the earlier, marked superseded, enters through its
        # successor alone, and still offers its own curve.
        write(
            joinpath(experiments, "retrieval.toml"),
            record(
                "correlated_with = [\"20000002\"]\ncorrelation_relation = \"republication\"";
                qualifiers = "[\"superseded: by 20000002\"]",
                second = "correlated_with = [\"20000001\"]\ncorrelation_relation = \"republication\"",
            ),
        )
        republished = run_pipeline(load_configuration(path; data_directory = directory))
        @test superseded_datasets(republished) == Dict("A. First 1979" => "B. First 1979")
        @test isempty(superseded_datasets(together))
        @test pooled_datasets(republished) == ["B. First 1979", "C. Other 2000"]
        @test isempty(republished.correlation_groups)
        @test republished.consensus_r_ν.ratio ≈ consensus([second, other]).ratio
        @test count(c -> c.kind == "dataset", republished.segmented_curves) == 3
        @test [only(entry.datasets) for entry in republished.leave_one_out] == ["B. First 1979", "C. Other 2000"]
        table = CSV.read(
            write_results(republished, joinpath(directory, "output", "republished"))["dataset_diagnostics"],
            DataFrame,
        )
        earlier = only(filter(r -> r.dataset == "A. First 1979", eachrow(table)))
        @test earlier.pooled == false
        @test earlier.exclusion_reason == "superseded by B. First 1979"
        @test run_metadata(republished)["result"]["superseded_datasets"] ==
            Dict("A. First 1979" => "B. First 1979")
    end
    return nothing
end

function test_flagged_curves(directory, configuration, result, written)
    @testset "a curve without a minimum over a range covering a window is flagged" begin
        # `partial` begins at A_H = 131, above the minimum at 130: its curve only rises. Without
        # a required window nothing says where a minimum lies, and nothing is flagged.
        partial = only(filter(c -> c.label == "partial", result.segmented_curves))
        @test occursin("rises from A_H = 131", unresolved_minimum(partial.fit))
        @test isempty(result.curve_flags)
        # A window the curve begins above: limited in range, which is not a flag.
        @test unresolved_minimum(partial.fit, [128:132]) === nothing
        @test occursin("covers the window 132:134", unresolved_minimum(partial.fit, [132:134]))
        full = only(filter(c -> c.label == "full", result.segmented_curves))
        @test unresolved_minimum(full.fit, [128:132]) === nothing

        # With the window required, the run flags the curve and still offers it.
        path = joinpath(directory, "windowed.toml")
        write(
            path,
            replace(
                PIPELINE_CONFIGURATION,
                "max_segments = 3" => "max_segments = 3\nrequired_windows = [[132, 134]]",
            ),
        )
        windowed = run_pipeline(load_configuration(path; data_directory = directory))
        @test collect(keys(windowed.curve_flags)) == ["partial"]
        reason = windowed.curve_flags["partial"]
        @test occursin("rises from A_H = 131", reason)
        @test occursin("covers the window 132:134", reason)
        @test any(c -> c.label == "partial", windowed.segmented_curves)
        written_ = write_results(windowed, joinpath(directory, "output", "windowed"))
        table = CSV.read(written_["dataset_diagnostics"], DataFrame)
        row(label) = only(filter(r -> r.dataset == label, eachrow(table)))
        @test row("partial").curve_flagged
        @test row("partial").curve_flag_reason == reason
        @test !row("full").curve_flagged
        @test !row("sparse").curve_flagged
        @test run_metadata(windowed)["result"]["flagged_curves"] == Dict("partial" => reason)
        @test isempty(run_metadata(result)["result"]["flagged_curves"])
        manifest = read_temperature_ratio_manifest(written_["manifest"])
        @test "partial" in curve_labels(manifest)
    end
    return nothing
end

@testset "pipeline" begin
    mktempdir() do directory
        # Vanishing shell corrections, with eq. (20) throughout, make a ∝ A, so every
        # fragmentation of A_H has a_L/a_H = (A₀ - A_H)/A_H and R_T follows in closed form. The
        # reader drops trailing zeros of S(Z) as padding, so the last row, at a nucleon number no
        # fragment carries, is not zero.
        open(joinpath(directory, "shell_corrections.dat"), "w") do io
            println(io, "n S_N S_Z")
            for n in 20:119
                println(io, n, " 0.0 0.0")
            end
            println(io, 120, " 1.0 1.0")
        end
        mkpath(joinpath(directory, "datasets"))
        write_multiplicity(joinpath(directory, "datasets", "full.dat"), 126:140)
        write_multiplicity(joinpath(directory, "datasets", "partial.dat"), 131:140)
        # Four pairs over fifteen mass numbers: below the coverage floor of 0.3.
        write_multiplicity(joinpath(directory, "datasets", "sparse.dat"), 126:4:138)
        path = joinpath(directory, "configuration.toml")
        write(path, PIPELINE_CONFIGURATION)

        configuration = load_configuration(path; data_directory = directory)
        result = run_pipeline(configuration)
        run_directory = joinpath(directory, "output", "run")
        written = write_results(result, run_directory)
        curve(label) = only(filter(c -> c.label == label, result.segmented_curves))

        @testset "the pin is applied at the symmetric split and nowhere else" begin
            full = curve("full")
            @test full.fit.x₀ == 126
            @test full.fit.pinned_value == 0.5
            @test first(full.r_ν.ratio) == 0.5
            @test first(full.R_T.ratio) ≈ 1 atol = 1e-12

            # The first complete pair of this dataset is A_H = 131, where r_ν ≈ 0.39. Pinning it
            # there would assert one half at a mass number where no identity holds.
            partial = curve("partial")
            @test partial.fit.x₀ == 131
            @test partial.fit.pinned_value === nothing
            @test first(partial.r_ν.ratio) ≈ reference_ratio(131) atol = 0.01

            @test curve(SYSTEMATIC_TREND_LABEL).fit.pinned_value == 0.5
        end

        @testset "the diagnostics and the curve table state which curves carry the pin" begin
            @test sort(result.symmetry.pinned_curves) == sort(["full", SYSTEMATIC_TREND_LABEL])
            @test isempty(result.symmetry.warnings)
            @test result.symmetry.charge_set_invariant === true
            @test result.symmetry.R_T_at_symmetric_split ≈ 1 atol = 1e-12

            table = CSV.read(written["segmented_curves"], DataFrame)
            @test basename(written["segmented_curves"]) ==
                "segmented_curves_$(run_identifier(configuration)).csv"
            flags = Dict(
                String(row.label) => row.pinned_at_symmetric_split for row in eachrow(table)
            )
            @test flags ==
                Dict("full" => true, "partial" => false, SYSTEMATIC_TREND_LABEL => true)
        end

        @testset "the manifest reads back with the domain the run used" begin
            manifest = read_temperature_ratio_manifest(written["manifest"])
            @test basename(written["manifest"]) ==
                "manifest_$(run_identifier(configuration)).toml"
            @test manifest.domain == manifest_domain(result)
            @test manifest.domain == ManifestDomain(
                result.model, result.averaging, result.domain, result.charge, result.masses
            )
            @test manifest.domain.level_density_model == "GC"
            @test !manifest.domain.deformed_branch
            @test manifest.domain.ratio_averaging == "charge_resolved"
            @test !manifest.domain.excitation_weighted
            @test manifest.domain.charges_per_mass == 5
            @test startswith(manifest.domain.charge_model, "Wahl1988(CF252S")
            @test manifest.domain.mass_table == "mass_excess_ame2020.dat"
            # One row of the curve table per manifest curve, keyed by its label.
            table = CSV.read(written["segmented_curves"], DataFrame)
            @test String.(table.label) == curve_labels(manifest)
            @test SYSTEMATIC_TREND_LABEL in curve_labels(manifest)
            @test manifest_curve(manifest, SYSTEMATIC_TREND_LABEL).kind == "systematic_trend"
            # The manifest holds the shared writer's fields and nothing besides; a tabulation
            # from no archive records no accession.
            document = TOML.parsefile(written["manifest"])
            @test Set(keys(document)) == Set(["system", "run", "domain", "segmented_curve"])
            for entry in document["segmented_curve"]
                @test Set(keys(entry)) == Set([
                    "label", "kind", "temperature_ratio_file", "multiplicity_ratio_pivots_file"
                ])
            end
        end

        @testset "a dataset below the coverage floor is pooled but offers no curve" begin
            @test !any(c -> c.label == "sparse", result.segmented_curves)
            @test startswith(result.dataset_outcomes["sparse"], "coverage")
            @test result.dataset_outcomes["full"] == "segmented curve"
            @test only(filter(d -> d.label == "sparse", result.dataset_diagnostics)).coverage ≈
                4 / 15
            # Pooled: its mass numbers reach the combined curve.
            @test 126 in result.consensus_r_ν.A_H
            # Diagnosed: the outcome is written, so the absence of a curve is explained.
            table = CSV.read(written["dataset_diagnostics"], DataFrame)
            @test nrow(table) == 3
            row = only(filter(r -> r.dataset == "sparse", eachrow(table)))
            @test startswith(row.segmented_curve, "coverage")
            @test row.coverage ≈ 4 / 15 atol = 1e-5
        end

        @testset "an exclusion that names no dataset read is refused" begin
            absent = Configuration(
                (
                    if name === :excluded_datasets
                        Dict("99999999" => "x")
                    else
                        getfield(configuration, name)
                    end for name in fieldnames(Configuration)
                )...
            )
            # Refused before anything is fitted: the log of the run holds no curve yet.
            logger = Test.TestLogger(; min_level = Base.CoreLogging.Info)
            @test_throws "which names no dataset read" Base.CoreLogging.with_logger(logger) do
                return run_pipeline(absent)
            end
            messages = [string(record.message) for record in logger.logs]
            @test any(contains("reading input"), messages)
            @test !any(contains("segmented curve"), messages)
            # A tabulation from no archive is excluded under its label.
            by_label = Configuration(
                (
                    if name === :excluded_datasets
                        Dict("sparse" => "few pairs")
                    else
                        getfield(configuration, name)
                    end for name in fieldnames(Configuration)
                )...
            )
            run = run_pipeline(by_label)
            # Named apart from the outer `written`, which the testsets below read.
            excluding_written = write_results(run, joinpath(directory, "output", "excluding"))
            table = CSV.read(excluding_written["dataset_diagnostics"], DataFrame)
            row(label) = only(filter(r -> r.dataset == label, eachrow(table)))
            @test row("sparse").pooled == false
            @test row("sparse").exclusion_reason == "few pairs"
            @test row("full").pooled == true
            @test row("partial").pooled == true
            @test pooled_datasets(run) == ["full", "partial"]
            @test pooled_datasets(result) == ["full", "partial", "sparse"]
        end

        test_trend_covariance(directory, configuration, result, written)
        test_leave_one_out(directory, configuration, result, written)
        test_experiments_pooled_as_one(directory, configuration, result, written)
        test_flagged_curves(directory, configuration, result, written)

        @testset "a yield exclusion that names no distribution is refused at the start" begin
            yields = joinpath(directory, "yields-check")
            mkpath(yields)
            write(
                joinpath(yields, "10000001_A.Both_2000.dat"),
                "A Y Y_uncertainty\n126 1.0 0.1\n134 2.0 0.1\n",
            )
            path = joinpath(directory, "yields-check.toml")
            write(
                path,
                PIPELINE_CONFIGURATION *
                "\n[yield]\nsubdirectory = \"yields-check\"\nmass_yield_file = \"yields-check/10000001_A.Both_2000.dat\"\n",
            )
            loaded = load_configuration(path; data_directory = directory)
            absent = Configuration(
                (
                    if name === :excluded_mass_yields
                        Dict("99999999" => "x")
                    else
                        getfield(loaded, name)
                    end for name in fieldnames(Configuration)
                )...
            )
            logger = Test.TestLogger(; min_level = Base.CoreLogging.Info)
            @test_throws "which names no distribution" Base.CoreLogging.with_logger(logger) do
                return run_pipeline(absent)
            end
            @test !any(contains("segmented curve"), [string(r.message) for r in logger.logs])

            # A configuration built by hand with an exclusion and no directory to apply it to
            # would write to the run record an exclusion that was never made.
            inapplicable = Configuration(
                (
                    if name === :excluded_mass_yields
                        Dict("10000001" => "x")
                    elseif name === :yield_directory
                        nothing
                    else
                        getfield(loaded, name)
                    end for name in fieldnames(Configuration)
                )...
            )
            @test_throws "without yield.subdirectory" run_pipeline(inapplicable)
            # And one whose directory is not there is refused as an argument, not by the file
            # system.
            misplaced = Configuration(
                (
                    if name === :yield_directory
                        joinpath(directory, "no-such-directory")
                    else
                        getfield(loaded, name)
                    end for name in fieldnames(Configuration)
                )...,
            )
            @test_throws ArgumentError run_pipeline(misplaced)
        end

        @testset "a dataset that forms no pair is not pooled, excluded or not" begin
            unpaired = joinpath(directory, "unpaired")
            mkpath(unpaired)
            write_multiplicity(joinpath(unpaired, "full.dat"), 126:140)
            # Heavy masses alone, as a coarse grid gives them: no mass has its complement.
            open(joinpath(unpaired, "grid.dat"), "w") do io
                println(io, "A nu nu_uncertainty")
                for A in 129:4:137
                    println(io, A, " 1.5 0.05")
                end
            end
            path = joinpath(directory, "unpaired.toml")
            body = replace(
                PIPELINE_CONFIGURATION,
                "subdirectory = \"datasets\"" => "subdirectory = \"unpaired\"",
            )
            write(path, body)
            run = run_pipeline(load_configuration(path; data_directory = directory))
            @test pooled_datasets(run) == ["full"]
            @test startswith(run.dataset_outcomes["grid"], "no complete fragment pair")
            table = CSV.read(
                write_results(run, joinpath(directory, "output", "unpaired"))["dataset_diagnostics"],
                DataFrame,
            )
            row(label) = only(filter(r -> r.dataset == label, eachrow(table)))
            @test row("grid").pooled == false
            @test ismissing(row("grid").exclusion_reason)
            @test row("full").pooled == true

            # Named by an exclusion, it stays out for the stated reason, and nothing else moves.
            write(
                path,
                replace(
                    body,
                    "[multiplicity]" => "[multiplicity]\nexclude = [{ dataset = \"grid\", reason = \"no pair on its grid\" }]",
                ),
            )
            excluding = run_pipeline(load_configuration(path; data_directory = directory))
            @test pooled_datasets(excluding) == ["full"]
            @test systematic_trend(excluding).fit.coefficients ==
                systematic_trend(run).fit.coefficients
            table = CSV.read(
                write_results(excluding, joinpath(directory, "output", "unpaired-excluded"))["dataset_diagnostics"],
                DataFrame,
            )
            @test row("grid").pooled == false
            @test row("grid").exclusion_reason == "no pair on its grid"
        end

        @testset "the run directory is written once and never into" begin
            @test isfile(joinpath(run_directory, "dataset_diagnostics.csv"))
            @test isfile(joinpath(run_directory, "r_nu_vs_A_H_pivots_full.csv"))
            # A tabulation from no archive has no accession, and its files keep its name.
            @test all(isempty, values(curve_accessions(result)))
            @test !any(contains("Cf252"), readdir(run_directory))
            @test_throws ArgumentError write_results(result, run_directory)
            # A run without yields reports the range mean and no total average.
            @test isempty(result.total_average_R_T)
            @test !haskey(written, "total_average_R_T")
            @test isfile(written["segmented_curves"])
            metadata = run_metadata(result)
            @test metadata["result"]["dataset_outcomes"]["sparse"] ==
                result.dataset_outcomes["sparse"]
            @test metadata["identifier"]["tokens"]["cov"] == 0.3
            @test !haskey(metadata["identifier"]["tokens"], "Ycov")
            @test metadata["configuration"]["min_pair_coverage"] == 0.3
            @test metadata["configuration"]["min_package_version"] == "0.2.3"
            # The version of the code that ran, as its project file states it now.
            @test metadata["source"]["package_version"] ==
                string(pkgversion(FissionTemperatureRatio))
            @test metadata["inputs"]["multiplicity_directory"] == "datasets"
        end

        @testset "R_T follows the closed form of the synthetic level density ratio" begin
            full = curve("full")
            for (index, A_H) in enumerate(full.R_T.A_H)
                r = full.r_ν.ratio[index]
                @test full.R_T.ratio[index] ≈ sqrt((1 - r) / (r * (252 - A_H) / A_H)) rtol =
                    1e-12
            end
        end

        @testset "yield distributions: one wing stands for the other, partial ones are not averaged" begin
            yields = joinpath(directory, "yields")
            mkpath(yields)
            peak(A_H) = exp(-((A_H - 134) / 4)^2)
            function write_yield(name, masses)
                open(joinpath(yields, name), "w") do io
                    println(io, "A Y Y_uncertainty")
                    for A in masses
                        println(io, A, " ", peak(max(A, 252 - A)), " 0.001")
                    end
                end
            end
            write_yield("10000001_A.Both_2000.dat", vcat(112:126, 127:140))
            write_yield("10000002_B.Light_2000.dat", 112:126)
            # Three heavy masses of fifteen: below the coverage floor of 0.3.
            write_yield("10000003_C.Tail_2000.dat", 138:140)
            path = joinpath(directory, "with_yields.toml")
            write(
                path,
                PIPELINE_CONFIGURATION *
                "\n[yield]\nsubdirectory = \"yields\"\nmass_yield_file = \"yields/10000001_A.Both_2000.dat\"\n",
            )
            with_yields = run_pipeline(load_configuration(path; data_directory = directory))

            @test with_yields.mass_yield_coverage["A. Both 2000"] == 1.0
            @test with_yields.mass_yield_coverage["B. Light 2000"] == 1.0
            # Coverage in the reference's yield: three tail masses hold a tenth of it.
            @test with_yields.mass_yield_coverage["C. Tail 2000"] ≈
                sum(peak, 138:140) / sum(peak, 126:140)
            trend = with_yields.total_average_R_T[SYSTEMATIC_TREND_LABEL]
            @test !haskey(trend, "C. Tail 2000")
            # The light wing measured alone gives the same heavy-fragment yields as both wings.
            @test trend["B. Light 2000"].value ≈ trend["A. Both 2000"].value rtol = 1e-12
            @test trend["A. Both 2000"].yield_fraction ≈ 1
            # A curve whose pairs begin at 131 takes in only the yield from there.
            partial = with_yields.total_average_R_T["partial"]["A. Both 2000"]
            @test partial.yield_fraction ≈ sum(peak, 131:140) / sum(peak, 126:140)
            metadata = run_metadata(with_yields)
            @test metadata["identifier"]["tokens"]["Ycov"] == 0.3
            @test collect(keys(metadata["result"]["mass_yields_not_averaged"])) ==
                ["C. Tail 2000"]
            @test startswith(
                metadata["result"]["mass_yields_not_averaged"]["C. Tail 2000"], "coverage"
            )
            table = CSV.read(
                write_results(with_yields, joinpath(directory, "output", "yields"))["total_average_R_T"],
                DataFrame,
            )
            @test !("C. Tail 2000" in table.mass_yield)
            @test all(0 .< table.yield_fraction .<= 1)

            # A distribution excluded by configuration is read and reported, never averaged.
            write(
                path,
                PIPELINE_CONFIGURATION *
                "\n[yield]\nsubdirectory = \"yields\"\nmass_yield_file = \"yields/10000001_A.Both_2000.dat\"\nexclude = [{ accession = \"10000002\", reason = \"not inclusive\" }]\n",
            )
            excluding = run_pipeline(load_configuration(path; data_directory = directory))
            @test haskey(excluding.mass_yield_coverage, "B. Light 2000")
            @test !any(
                haskey(per, "B. Light 2000") for per in values(excluding.total_average_R_T)
            )
            @test run_metadata(excluding)["result"]["mass_yields_not_averaged"]["B. Light 2000"] ==
                "excluded by configuration: not inclusive"
            @test run_metadata(excluding)["identifier"]["hashed"]["Y"]["exclude"] ==
                ["10000002"]

            # The two floors are separate keys: lowering one leaves the other gate where it was.
            write(
                path,
                PIPELINE_CONFIGURATION *
                "\n[yield]\nsubdirectory = \"yields\"\nmass_yield_file = \"yields/10000001_A.Both_2000.dat\"\nmin_yield_coverage = 0.05\n",
            )
            lower_yield = run_pipeline(load_configuration(path; data_directory = directory))
            @test haskey(lower_yield.total_average_R_T[SYSTEMATIC_TREND_LABEL], "C. Tail 2000")
            @test startswith(lower_yield.dataset_outcomes["sparse"], "coverage")
            write(
                path,
                replace(
                    PIPELINE_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nmin_pair_coverage = 0.2",
                ) *
                "\n[yield]\nsubdirectory = \"yields\"\nmass_yield_file = \"yields/10000001_A.Both_2000.dat\"\n",
            )
            lower_pair = run_pipeline(load_configuration(path; data_directory = directory))
            @test !haskey(lower_pair.total_average_R_T[SYSTEMATIC_TREND_LABEL], "C. Tail 2000")
            @test !startswith(lower_pair.dataset_outcomes["sparse"], "coverage")
        end

        @testset "the trend of a pool of one interpolated dataset is that dataset's fit" begin
            single = joinpath(directory, "single")
            mkpath(single)
            file = "10000009_A.Interp_2001.dat"
            write_multiplicity(joinpath(single, file), 126:140)
            # 29 rows written from 20 measured masses: every point counts for 20/29 of one.
            write(
                joinpath(single, "retrieval.toml"),
                """
                [[accepted]]
                file = "$(file)"
                identifier = "10000009"
                mass_treatment = "interpolated"
                mass_values_non_integer = 20
                rows_written = 29
                qualifiers = []

                [run]
                package_version = "0.2.3"
                """,
            )
            path = joinpath(directory, "single.toml")
            write(
                path,
                replace(
                    PIPELINE_CONFIGURATION,
                    "subdirectory = \"datasets\"" => "subdirectory = \"single\"",
                ),
            )
            one = run_pipeline(load_configuration(path; data_directory = directory))
            dataset = only(filter(c -> c.kind == "dataset", one.segmented_curves)).fit
            trend = systematic_trend(one).fit
            @test dataset.measured_points ≈ 15 * 20 / 29
            @test segments(trend) > 1
            # The pooling fraction enters once: the same weights, χ², degrees of freedom,
            # selection and covariance as the dataset's own fit, to the last bit.
            for field in
                (:breakpoints, :coefficients, :wrss, :dof, :bic, :selection, :covariance)
                @test getfield(trend, field) == getfield(dataset, field)
            end
            # The label is the display name; the accession is in the record and in the names of
            # the files.
            @test curve_accessions(one) ==
                Dict("A. Interp 2001" => "10000009", SYSTEMATIC_TREND_LABEL => "")
            files = readdir(
                dirname(write_results(one, joinpath(directory, "output", "one"))["manifest"])
            )
            for name in ("R_T_vs_A_H_segmented", "r_nu_vs_A_H_pivots", "r_nu_vs_A_H")
                @test "$(name)_10000009_A.Interp_2001.csv" in files
            end
            @test "R_T_vs_A_H_segmented_systematic_trend.csv" in files
            manifest = read_temperature_ratio_manifest(
                staged_manifest(joinpath(directory, "output", "one"))
            )
            @test manifest_curve(manifest, "A. Interp 2001").accession == "10000009"
            @test manifest_curve(manifest, SYSTEMATIC_TREND_LABEL).accession == ""
            # The combined curve states the standard error, which carries the fraction.
            @test all(one.consensus_r_ν.σ .≈ only(one.r_ν).σ ./ sqrt(20 / 29))
        end
    end
end
