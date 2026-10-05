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
# Both need the measured input, so both are reported as skipped on a bare clone. Where the input
# of a system is held, every row of that system must be found: a dataset or a distribution the
# archive no longer returns takes its row out of `fixtures.jl` by an edit that states why, not by
# a skip.

function published_rows_hold(results; widened = Dict())
    for (system, dataset, distribution, published, tolerance) in PUBLISHED_TOTAL_AVERAGES
        # The input of this system is not held: nothing to compare.
        if !haskey(results, system)
            @test_skip false
            continue
        end
        tolerance = get(widened, (system, dataset, distribution), tolerance)
        averages = results[system].total_average_R_T
        held = haskey(averages, dataset) && haskey(averages[dataset], distribution)
        held || @error "no total average for a row of Table 1" system dataset distribution
        @test held
        held && @test averages[dataset][distribution].value ≈ published rtol = tolerance
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
                    results[system] = run_pipeline(configuration; leave_one_out = false)
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
                # The refits of the trend with a dataset left out, for the smallest pool only.
                results[system] = run_pipeline(
                    configuration; leave_one_out = system == "U233_nth"
                )
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
            ("Pu239_nth", "23012008", true),
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

        # The covariance of every trend carries the autocorrelation of the pooled datasets'
        # deviations, and its scale refers to what a fit to such errors leaves in the residuals.
        for run in values(results)
            trend = systematic_trend(run)
            @test 0.6 < deviation_autocorrelation(run) < 0.95
            @test trend.fit.correlated
            @test trend.fit.expected_wrss < trend.fit.dof
            @test all(
                !c.fit.correlated && c.fit.expected_wrss == c.fit.dof for
                c in run.segmented_curves if c.kind == "dataset"
            )
        end

        # A dataset curve that resolves no minimum is flagged and stays offered; a large χ² under
        # small quoted uncertainties flags nothing.
        if haskey(results, "Pu239_nth")
            run = results["Pu239_nth"]
            @test occursin("no interior minimum", run.curve_flags["C. Tsuchiya 2000"])
            @test any(c -> c.label == "C. Tsuchiya 2000", run.segmented_curves)
        end
        if haskey(results, "Cf252_sf")
            # Bowman, Britt and Mehta begin above the window: limited in range, not flagged.
            @test isempty(results["Cf252_sf"].curve_flags)
        end
        if haskey(results, "U233_nth")
            @test collect(keys(results["U233_nth"].curve_flags)) == ["J.S. Fraser 1966"]
        end
        haskey(results, "U235_nth") && @test isempty(results["U235_nth"].curve_flags)

        # Basova 1979 and Zamyatnin 1979 are two reductions of one experiment: pooled as one.
        for system in ("Cf252_sf", "Pu239_nth")
            haskey(results, system) || continue
            @test results[system].correlation_groups ==
                [["Yu.S. Zamyatnin 1979", "B.G. Basova 1979"]]
        end

        # With one pooled dataset left out in turn: one refit per pooled dataset, and the
        # jackknife uncertainty of the trend's ⟨R_T⟩ over the primary distribution.
        if haskey(results, "U233_nth")
            run = results["U233_nth"]
            @test reduce(vcat, [entry.datasets for entry in run.leave_one_out]) == pooled_datasets(run)
            @test all(entry.outcome == "segmented curve" for entry in run.leave_one_out)
            spread = leave_one_out_spread(run, "V.M. Surin 1972")
            trend = run.total_average_R_T[SYSTEMATIC_TREND_LABEL]["V.M. Surin 1972"]
            @test spread.min < trend.value < spread.max
            @test spread.uncertainty > 0
        end
        haskey(results, "Cf252_sf") && @test isempty(results["Cf252_sf"].leave_one_out)

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
