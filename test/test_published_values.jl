# Regression against the total averages ⟨R_T⟩ of Table 1 of Eur. Phys. J. A 60, 190 (2024),
# doi:10.1140/epja/s10050-024-01375-7; the rows and their tolerances are in `fixtures.jl`.
#
# Two runs. The published settings — the ratio of means over the evaluated Gaussian charge tables
# with ΔZ(A₀/2) = 0, or over the mean values for 233-U — need the digitised tables, which are held
# locally and not shipped, so that run is skipped without them. The shipped settings — the
# charge-resolved inversion on Wahl's 1988 model, weighted by the ⟨TXE⟩ of the configured
# ⟨TKE⟩(A) dataset, and Y(A) symmetrized to the pre-neutron identity — are held to the same
# tolerances, seven of them widened (`SHIPPED_TOLERANCE`). Both runs average over every
# distribution held, with the primary experiment's as the coverage reference: the shipped
# configurations average over the primary alone, which for 233-U, 235-U and 239-Pu is not one the
# table names. The 239-Pu rows are compared over Nishio 1995, a yield distribution inferred by
# agreement, since the table's caption names none.
#
# Both need the measured input, so both are reported as skipped on a bare clone. A row is skipped, not failed,
# when the local data holds no dataset of that label: what the archive returns for a query changes
# over time.

function published_rows_hold(results; widened = Dict())
    for (system, dataset, distribution, published, tolerance) in PUBLISHED_TOTAL_AVERAGES
        tolerance = get(widened, (system, dataset, distribution), tolerance)
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

# Without the measured input, which is not shipped, the testset is reported as skipped.
DATA_AVAILABLE || @testset "published total averages" begin
    @test_skip DATA_AVAILABLE
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

    @testset "with the shipped settings" begin
        results = Dict{String,Any}()
        mktempdir() do directory
            for system in filter(has_inputs, PUBLISHED_SYSTEMS)
                configuration = variant_configuration(
                    system,
                    directory;
                    set = Dict("yield" => Dict("subdirectory" => "$(system)/Y_vs_A")),
                )
                results[system] = run_pipeline(configuration)
                domain = manifest_domain(results[system])
                @test domain.ratio_averaging == "charge_resolved"
                @test domain.excitation_weighted
                @test startswith(domain.charge_model, "Wahl1988(")
            end
        end
        published_rows_hold(results; widened = SHIPPED_TOLERANCE)

        # The shipped exclusions, by accession, each with its reason on record. Zeynalov 2019 is
        # kept out of the 252-Cf pool (41739002) and of the 235-U pool (41738002) and still
        # offers its own curve; Batenkov 2004 (41502005 for 235-U, 41502006 for 239-Pu) forms no
        # fragment pair on its 4-u grid and offers none. The datasets must be held: an exclusion
        # of an absent dataset is refused when the configuration is loaded, and `only` fails
        # where no dataset has the accession.
        for (system, accession, offers_curve) in (
            ("Cf252_sf", "41739002", true),
            ("U235_nth", "41738002", true),
            ("U235_nth", "41502005", false),
            ("Pu239_nth", "41502006", false),
        )
            haskey(results, system) || continue
            run = results[system]
            label = only(l for (l, a) in curve_accessions(run) if a == accession)
            @test run.configuration.excluded_datasets[accession] isa String
            @test !(label in pooled_datasets(run))
            @test any(c -> c.label == label, run.segmented_curves) == offers_curve
            if offers_curve
                @test run.dataset_outcomes[label] == "segmented curve"
            else
                @test startswith(run.dataset_outcomes[label], "no complete fragment pair")
            end
            excluded = run_metadata(run)["configuration"]["excluded_datasets"]
            @test occursin(accession, excluded[accession])
        end

        # Straede's 235-U yields are spectrum-averaged: used, and flagged. Only a distribution the
        # run read is reported, whatever else the retrieval record lists.
        if haskey(results, "U235_nth")
            flagged = run_metadata(results["U235_nth"])["inputs"]["flagged_mass_yields"]
            @test flagged["Ch.Straede 1987"] == ["SPA"]
            read = Set(y.label for y in results["U235_nth"].mass_yields)
            @test issubset(keys(flagged), read)
        end
    end
end
