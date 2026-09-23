# Regression against the total averages ⟨R_T⟩ of Table 1 of Eur. Phys. J. A 60, 190 (2024),
# doi:10.1140/epja/s10050-024-01375-7, back-shifted Fermi gas, five charges per mass number.
#
# Needs the input data and the shipped configurations, so it is skipped on a bare clone. A row is
# skipped, not failed, when the local data holds no dataset of that label: what the archive returns
# for a query changes over time.
#
# Tolerances are relative and per row. The inputs were retrieved independently of whatever the
# authors used, and the breakpoints are selected rather than placed by hand, so agreement is at the
# per-cent level, not to the digits quoted. Fraser (233-U) quotes no uncertainties and has fifteen
# usable pairs; the published value itself carries ±0.25.

const PUBLISHED_TOTAL_AVERAGES = [
    # system, ν(A) dataset, Y(A) distribution, published ⟨R_T⟩, relative tolerance
    ("U233_nth", "K. Nishio 1998", "V.M. Surin 1972", 1.1861, 0.006),
    ("U233_nth", "V.F. Apalin 1965", "V.M. Surin 1972", 1.0346, 0.008),
    ("U233_nth", "J.S. Fraser 1966", "V.M. Surin 1972", 1.3809, 0.015),
    ("Cf252_sf", "C. Budtz-jorgensen 1988", "A. Goeoek 2014", 1.0975, 0.006),
    ("Cf252_sf", "Yu.S. Zamyatnin 1979", "A. Goeoek 2014", 1.1128, 0.006),
    ("Cf252_sf", "A. Goeoek 2014", "A. Goeoek 2014", 1.1163, 0.008),
    ("Cf252_sf", "A. Al-adili 2020", "A. Goeoek 2014", 1.0860, 0.008),
    ("U235_nth", "K. Nishio 1998", "A. Al-adili 2020", 1.1586, 0.008),
    ("U235_nth", "K. Nishio 1998", "Ch.Straede 1987", 1.1644, 0.012),
    ("U235_nth", "A.S. Vorobyev 2010", "A. Al-adili 2020", 1.1186, 0.010),
    ("U235_nth", "A.S. Vorobyev 2010", "Ch.Straede 1987", 1.1221, 0.006),
]

DATA_AVAILABLE && @testset "published total averages" begin
    systems = unique(first.(PUBLISHED_TOTAL_AVERAGES))
    results = Dict{String,Any}()
    for system in systems
        all(isdir(joinpath(DATA_DIRECTORY, system, d)) for d in ("nu_vs_A", "Y_vs_A")) ||
            continue
        configuration = load_configuration(
            joinpath(pkgdir(FissionTemperatureRatio), "config", "$(system).toml");
            data_directory = DATA_DIRECTORY,
        )
        results[system] = run_pipeline(configuration)
    end

    for (system, dataset, distribution, published, tolerance) in PUBLISHED_TOTAL_AVERAGES
        averages = haskey(results, system) ? results[system].total_average_R_T : Dict()
        if haskey(averages, dataset) && haskey(averages[dataset], distribution)
            @test averages[dataset][distribution].value ≈ published rtol = tolerance
        else
            @test_skip false
        end
    end
end
