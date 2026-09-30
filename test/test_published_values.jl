# Regression against the total averages ⟨R_T⟩ of Table 1 of Eur. Phys. J. A 60, 190 (2024),
# doi:10.1140/epja/s10050-024-01375-7; the rows and their tolerances are in `fixtures.jl`.
#
# Two runs. The published settings — the ratio of means over the evaluated Gaussian charge tables
# with ΔZ(A₀/2) = 0, or over the mean values for 233-U — need the digitised tables, which are held
# locally and not shipped, so that run is skipped without them. The shipped configurations — the
# charge-resolved inversion on Wahl's 1988 model, weighted by ⟨TXE⟩ where a ⟨TKE⟩(A) dataset is
# configured — are held to the same tolerances.
#
# Both need the measured input, so both are skipped on a bare clone. A row is skipped, not failed,
# when the local data holds no dataset of that label: what the archive returns for a query changes
# over time.

function published_rows_hold(results)
    for (system, dataset, distribution, published, tolerance) in PUBLISHED_TOTAL_AVERAGES
        averages = haskey(results, system) ? results[system].total_average_R_T : Dict()
        if haskey(averages, dataset) && haskey(averages[dataset], distribution)
            @test averages[dataset][distribution].value ≈ published rtol = tolerance
        else
            @test_skip false
        end
    end
end

const PUBLISHED_SYSTEMS = unique(first.(PUBLISHED_TOTAL_AVERAGES))

function has_inputs(system)
    return all(isdir(joinpath(DATA_DIRECTORY, system, d)) for d in ("nu_vs_A", "Y_vs_A"))
end

DATA_AVAILABLE && @testset "published total averages" begin
    @testset "with the published settings" begin
        if PUBLISHED_TABLES_AVAILABLE
            results = Dict{String,Any}()
            mktempdir() do directory
                for system in filter(has_inputs, PUBLISHED_SYSTEMS)
                    settings = published_settings(system)
                    configuration = variant_configuration(
                        system, directory; set = settings.set, remove = settings.remove
                    )
                    results[system] = run_pipeline(configuration)
                    @test manifest_domain(results[system]).ratio_averaging == "ratio_of_means"
                end
            end
            published_rows_hold(results)
        else
            @test_skip false
        end
    end

    @testset "with the shipped configurations" begin
        results = Dict{String,Any}()
        for system in filter(has_inputs, PUBLISHED_SYSTEMS)
            configuration = load_configuration(
                joinpath(CONFIG_DIRECTORY, "$(system).toml");
                data_directory = DATA_DIRECTORY,
            )
            results[system] = run_pipeline(configuration)
            @test manifest_domain(results[system]).ratio_averaging == "charge_resolved"
            @test startswith(manifest_domain(results[system]).charge_model, "Wahl1988(")
        end
        published_rows_hold(results)

        # Straede's 235-U yields are spectrum-averaged: used, and flagged. The retrieval record
        # lists many more such sets; only a distribution the run read is reported.
        if haskey(results, "U235_nth")
            flagged = run_metadata(results["U235_nth"])["inputs"]["flagged_mass_yields"]
            @test flagged == Dict{String,Any}("Ch.Straede 1987" => ["SPA"])
        end
    end
end
