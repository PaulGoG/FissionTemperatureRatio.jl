# The handoff to a consuming code, exercised rather than assumed.
#
# Everything below reads only what a run writes, through the reader a consuming code uses:
# FissionFragmentsDomain's `staged_manifest`, `read_temperature_ratio_manifest`,
# `temperature_ratio_path` and `read_segmented_curve`. That is the contract a prompt emission code
# depends on; a change that breaks it would otherwise be found by whoever is downstream rather than
# here.

DATA_AVAILABLE && @testset "the handoff to a consuming code" begin
    configuration = load_configuration(
        joinpath(CONFIG_DIRECTORY, "U233_nth.toml"); data_directory = DATA_DIRECTORY
    )

    mktempdir() do root
        result = run_pipeline(configuration)
        # Staged as the consumer stages it: the whole run directory, copied.
        directory = joinpath(root, "sims", "U233_nth", run_identifier(configuration))
        write_results(result, directory)
        path = staged_manifest(directory)
        manifest = read_temperature_ratio_manifest(path)

        @testset "the consumer's reader accepts the manifest" begin
            @test basename(path) == "manifest_$(run_identifier(configuration)).toml"
            @test manifest.ordinate == "R_T"
            @test manifest.abscissa == ["A_H"]
            @test manifest.columns == ["A_H", "R_T", "R_T_uncertainty"]
            labels = curve_labels(manifest)
            @test !isempty(labels)
            @test count(==(SYSTEMATIC_TREND_LABEL), labels) == 1
            @test all(c -> c.kind in MANIFEST_CURVE_KINDS, manifest.curves)
        end

        @testset "the domain the run used is on record" begin
            @test manifest.domain == manifest_domain(result)
            @test manifest.domain.ratio_averaging == "charge_resolved"
            # The 233-U configuration names Geltenbort's ⟨TKE⟩(A), so the inversion is weighted.
            @test manifest.domain.excitation_weighted
            @test manifest.domain.level_density_model == "BSFG"
            @test !manifest.domain.deformed_branch
            @test manifest.domain.charge_model == charge_model_label(result.charge)
            @test startswith(manifest.domain.charge_model, "Wahl1988(U233T")
        end

        @testset "the system is identified without parsing a label" begin
            system = manifest.system
            @test system.label == "U233_nth"
            @test system.notation == "²³³U(nth,f)"
            @test (system.target_A, system.target_Z) == (233, 92)
            @test system.channel == "nth"
            @test system.reaction == "n,f"
            @test (system.compound_A, system.compound_Z) == (234, 92)
        end

        @testset "every curve resolves to a dense tabulation of R_T" begin
            for label in curve_labels(manifest)
                curve = read_segmented_curve(temperature_ratio_path(manifest, label))
                @test curve.masses == collect(first(curve.masses):last(curve.masses))
                @test all(>(0), curve.values)
                # Tabulated at every mass number, so the consumer's interpolation reproduces the
                # tabulated value exactly, and the curve is undefined outside its range.
                for (A, value) in zip(curve.masses, curve.values)
                    @test curve(A) == value
                end
                @test curve(first(curve.masses) - 1) === nothing
                @test curve(last(curve.masses) + 1) === nothing
            end
        end

        @testset "the exact identity at the symmetric split survives the round trip" begin
            A_symmetric = configuration.system.compound.A ÷ 2
            for label in curve_labels(manifest)
                curve = read_segmented_curve(temperature_ratio_path(manifest, label))
                first(curve.masses) == A_symmetric || continue
                @test curve(A_symmetric) ≈ 1 atol = 1e-10
            end
        end

        @testset "what the manifest points to is what the run produced" begin
            curve = read_segmented_curve(temperature_ratio_path(manifest, SYSTEMATIC_TREND_LABEL))
            produced = systematic_trend(result).R_T
            @test curve.masses == produced.A_H
            # The file is rounded to the configured number of significant figures.
            tolerance = 5 * 10.0^(-configuration.output.significant_digits)
            @test all(isapprox.(curve.values, produced.ratio; rtol = tolerance))
        end

        @testset "the per-curve table and the totals carry the run token" begin
            identifier = run_identifier(configuration)
            table = CSV.read(
                joinpath(directory, "segmented_curves_$(identifier).csv"), DataFrame
            )
            @test String.(table.label) == curve_labels(manifest)
            trend = only(filter(row -> row.label == SYSTEMATIC_TREND_LABEL, eachrow(table)))
            @test trend.kind == "systematic_trend"
            @test trend.first_A_H == first(systematic_trend(result).R_T.A_H)
            @test trend.pairs == length(result.consensus_r_ν)
            @test all(0 .< table.coverage .<= 1)

            averages = CSV.read(
                joinpath(directory, "total_average_R_T_$(identifier).csv"), DataFrame
            )
            @test names(averages) == [
                "segmented_curve",
                "mass_yield",
                "R_T",
                "R_T_uncertainty",
                "R_T_uncertainty_independent_points",
            ]
            @test all(averages.R_T_uncertainty .≥ 0)
            # The token is on the directory, the manifest and these two tables; the per-curve
            # ratio tables are named by quantity and dataset, and found through the manifest.
            @test count(contains(identifier), readdir(directory)) == 3
        end

        @testset "excluded, unfitted and qualified datasets are visible, not missing" begin
            table = CSV.read(joinpath(directory, "dataset_diagnostics.csv"), DataFrame)
            @test nrow(table) == length(result.datasets)
            @test all(in(("true", "false")), string.(table.pooled))
            @test all(in(("true", "false")), string.(table.flagged))
            @test Set(string.(table.dataset)) ⊇
                Set(c.label for c in manifest.curves if c.kind == "dataset")
            @test all(!isempty, string.(table.segmented_curve))
            files = readdir(directory)
            @test any(startswith("r_nu_vs_A_H_segmented_"), files)
            @test "r_nu_vs_A_H_consensus_systematic_trend.csv" in files
        end

        @testset "the run metadata records the ⟨TKE⟩(A) retrieval" begin
            metadata = run_metadata(result)
            record = metadata["mean_kinetic_energy"]
            @test record["excitation_weighted"]
            @test record["file"] == "U233_nth/TKE_vs_A/21981008_P.Geltenbort_1985.dat"
            @test haskey(record["retrieval_run"], "package_revision")
            @test length(record["retrieval_record_sha1"]) == 40
            # Geltenbort tabulates A_H to 158; the range reaches 159.
            @test record["extrapolated"] == [159]
            # Its yield-weighted mean against the standard, 170.1 MeV for 233-U(nth,f).
            for (_, offset) in record["offset_from_standard"]
                @test offset["standard_MeV"] == 170.1
                @test offset["offset_MeV"] ≈ offset["mean_MeV"] - 170.1
                @test abs(offset["offset_MeV"]) < 15
            end
            @test metadata["domain"]["excitation_weighted"]
        end
    end
end
