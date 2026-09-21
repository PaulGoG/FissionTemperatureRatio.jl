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
mass_excess_file = "mass_excess.dat"
shell_correction_file = "shell_corrections.dat"

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
        write(joinpath(directory, "mass_excess.dat"), "1 1 H 7289.0 0.0\n1 0 n 8071.0 0.0\n")
        # Vanishing shell corrections make a ∝ A, so R_a(A_H) = (A₀ - A_H)/A_H in closed form.
        open(joinpath(directory, "shell_corrections.dat"), "w") do io
            println(io, "n S_N S_Z")
            for n in 20:120
                println(io, n, " 0.0 0.0")
            end
        end
        mkpath(joinpath(directory, "datasets"))
        write_multiplicity(joinpath(directory, "datasets", "full.dat"), 126:140)
        write_multiplicity(joinpath(directory, "datasets", "partial.dat"), 131:140)
        path = joinpath(directory, "configuration.toml")
        write(path, PIPELINE_CONFIGURATION)

        configuration = load_configuration(path; data_directory = directory)
        result = run_pipeline(configuration; output_root = joinpath(directory, "output"))
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

            @test curve(TREND_LABEL).fit.pinned_value == 0.5
        end

        @testset "the diagnostics and the manifest state which curves carry the pin" begin
            @test sort(result.diagnostics["r_nu_pinned_at_symmetry"]) ==
                sort(["full", TREND_LABEL])
            @test !any(first(w) == "pin" for w in result.diagnostics["warnings"])

            results_directory = joinpath(directory, "output", "results", "Cf252_sf")
            manifest = TOML.parsefile(
                joinpath(
                    results_directory,
                    only(filter(startswith("manifest_"), readdir(results_directory))),
                ),
            )
            flags = Dict(
                entry["label"] => entry["pinned_at_symmetric_split"] for
                entry in manifest["segmented_curve"]
            )
            @test flags == Dict("full" => true, "partial" => false, TREND_LABEL => true)
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
