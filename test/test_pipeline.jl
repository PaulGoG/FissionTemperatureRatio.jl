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
            # The manifest holds the shared writer's fields and nothing besides.
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

        @testset "the run directory is written once and never into" begin
            @test isfile(joinpath(run_directory, "dataset_diagnostics.csv"))
            @test isfile(joinpath(run_directory, "r_nu_vs_A_H_pivots_full.csv"))
            @test !any(contains("Cf252"), readdir(run_directory))
            @test_throws ArgumentError write_results(result, run_directory)
            # A run without yields reports the range mean and no total average.
            @test isempty(result.total_average_R_T)
            @test !haskey(written, "total_average_R_T")
            @test isfile(written["segmented_curves"])
            metadata = run_metadata(result)
            @test metadata["result"]["dataset_outcomes"]["sparse"] ==
                result.dataset_outcomes["sparse"]
            @test metadata["identifier"]["tokens"]["mincov"] == 0.3
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
    end
end
