# The configuration parser must enforce exactly the constraints its comments document, and fail
# naming the offending key, so that a run cannot start from a configuration it cannot honour.

const MINIMAL_CONFIGURATION = """
[system]
target_A = 252
target_Z = 98
channel = "sf"

[fragmentation]
charges_per_mass = 5
heavy_mass_max = 140

[level_density]
model = "BSFG"

[multiplicity]
subdirectory = "datasets"

[segments]
max_segments = 3

[output]
significant_digits = 6
"""

"Write `body` as a second configuration beside `path`, and return its path."
function write_beside(path::AbstractString, body::AbstractString)
    other = joinpath(dirname(path), "other_$(length(readdir(dirname(path)))).toml")
    write(other, body)
    return other
end

function with_configuration(f, body::AbstractString)
    mktempdir() do directory
        write(joinpath(directory, "mass_excess.dat"), "1 1 H 7289.0 0.0\n1 0 n 8071.0 0.0\n")
        mkpath(joinpath(directory, "datasets"))
        # A ⟨TKE⟩(A) dataset with the record of the retrieval that wrote it, and one without.
        mkpath(joinpath(directory, "TKE_vs_A"))
        write(
            joinpath(directory, "TKE_vs_A", "1_A.Author_2000.dat"),
            "A TKE\n126 170.0\n140 175.0\n",
        )
        write(joinpath(directory, "TKE_vs_A", "2_B.Author_2001.dat"), "A TKE\n126 170.0\n")
        write(
            joinpath(directory, "TKE_vs_A", "retrieval.toml"),
            "[[accepted]]\nfile = \"1_A.Author_2000.dat\"\nqualifiers = []\n\n[run]\npackage_version = \"0.2.3\"\n",
        )
        path = joinpath(directory, "configuration.toml")
        write(path, body)
        return f(path, directory)
    end
end

@testset "configuration" begin
    @testset "a valid configuration parses" begin
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test configuration.system == spontaneous_fission(Nuclide(98, 252))
            @test configuration.system.compound.A == 252
            @test reaction(configuration.system) == "0,f"
            # The smallest heavy mass defaults to the symmetric split.
            @test configuration.fragmentation.heavy_mass_min == 126
            @test A_H_range(configuration) == 126:140
            @test has_symmetric_split(configuration)
            # Defaults: Wahl's charge model, the shipped mass table, the charge-resolved
            # inversion without excitation weights, both Gilbert-Cameron formulas.
            @test configuration.fragmentation.charge_distribution == "wahl"
            @test !configuration.fragmentation.zero_polarization_at_symmetry
            @test configuration.level_density.mass_excess_file == "ame2020"
            @test configuration.level_density.shell_correction_file == "gilbert_cameron_1965"
            @test configuration.level_density.ratio_averaging == "charge_resolved"
            @test configuration.level_density.mean_kinetic_energy_file === nothing
            @test configuration.level_density.deformed_branch
            @test configuration.symmetrize_yields
            @test configuration.segments.pin_symmetric_split
            @test configuration.segments.min_segment_span == 3
            @test configuration.segments.min_pair_coverage == 0.3
            @test configuration.min_yield_coverage == 0.3
            @test configuration.min_retrieval_version == v"0.2.3"
            @test configuration.output.significant_digits == 6
            @test configuration.data_directory == directory
        end
    end

    @testset "missing file" begin
        @test_throws ArgumentError load_configuration(joinpath(@__DIR__, "absent.toml"))
    end

    @testset "each constraint is enforced" begin
        cases = [
            (
                "missing section",
                replace(MINIMAL_CONFIGURATION, "[output]\nsignificant_digits = 6\n" => ""),
            ),
            ("missing key", replace(MINIMAL_CONFIGURATION, "target_A = 252\n" => "")),
            (
                "wrong type",
                replace(MINIMAL_CONFIGURATION, "target_A = 252" => "target_A = \"252\""),
            ),
            (
                "even charges per mass",
                replace(
                    MINIMAL_CONFIGURATION, "charges_per_mass = 5" => "charges_per_mass = 4"
                ),
            ),
            (
                "charge above mass",
                replace(MINIMAL_CONFIGURATION, "target_Z = 98" => "target_Z = 300"),
            ),
            (
                "unknown channel",
                replace(MINIMAL_CONFIGURATION, "channel = \"sf\"" => "channel = \"p,f\""),
            ),
            (
                "energy given for spontaneous fission",
                replace(
                    MINIMAL_CONFIGURATION,
                    "channel = \"sf\"" => "channel = \"sf\"\nincident_energy = 1.0",
                ),
            ),
            (
                "heavy mass range below symmetry",
                replace(
                    MINIMAL_CONFIGURATION, "heavy_mass_max = 140" => "heavy_mass_max = 100"
                ),
            ),
            (
                "heavy mass range above A0",
                replace(
                    MINIMAL_CONFIGURATION, "heavy_mass_max = 140" => "heavy_mass_max = 300"
                ),
            ),
            (
                "heavy_mass_min below the symmetric split",
                replace(
                    MINIMAL_CONFIGURATION,
                    "heavy_mass_max = 140" => "heavy_mass_min = 120\nheavy_mass_max = 140",
                ),
            ),
            (
                "heavy_mass_min above heavy_mass_max",
                replace(
                    MINIMAL_CONFIGURATION,
                    "heavy_mass_max = 140" => "heavy_mass_min = 141\nheavy_mass_max = 140",
                ),
            ),
            (
                "unknown level density model",
                replace(MINIMAL_CONFIGURATION, "model = \"BSFG\"" => "model = \"XYZ\""),
            ),
            (
                "absent shell correction table",
                replace(
                    MINIMAL_CONFIGURATION,
                    "model = \"BSFG\"" => "model = \"GC\"\nshell_correction_file = \"absent.dat\"",
                ),
            ),
            (
                "unknown averaging",
                replace(
                    MINIMAL_CONFIGURATION,
                    "model = \"BSFG\"" => "model = \"BSFG\"\nratio_averaging = \"other\"",
                ),
            ),
            (
                "absent mass excess file",
                replace(
                    MINIMAL_CONFIGURATION,
                    "model = \"BSFG\"" => "model = \"BSFG\"\nmass_excess_file = \"absent.dat\"",
                ),
            ),
            (
                "absent charge distribution table",
                replace(
                    MINIMAL_CONFIGURATION,
                    "heavy_mass_max = 140" => "heavy_mass_max = 140\ncharge_distribution_file = \"absent.dat\"",
                ),
            ),
            (
                "zero polarization imposed on Wahl's model",
                replace(
                    MINIMAL_CONFIGURATION,
                    "heavy_mass_max = 140" => "heavy_mass_max = 140\nzero_polarization_at_symmetry = true",
                ),
            ),
            (
                "Gilbert-Cameron branch not a Boolean",
                replace(
                    MINIMAL_CONFIGURATION,
                    "model = \"BSFG\"" => "model = \"BSFG\"\ndeformed_branch = 1",
                ),
            ),
            (
                "excitation weights on an effective ratio",
                replace(
                    MINIMAL_CONFIGURATION,
                    "model = \"BSFG\"" => "model = \"BSFG\"\nratio_averaging = \"ratio_of_means\"\nmean_kinetic_energy_file = \"TKE_vs_A/1_A.Author_2000.dat\"",
                ),
            ),
            (
                "absent kinetic energy file",
                replace(
                    MINIMAL_CONFIGURATION,
                    "model = \"BSFG\"" => "model = \"BSFG\"\nmean_kinetic_energy_file = \"TKE_vs_A/absent.dat\"",
                ),
            ),
            (
                "kinetic energy file no retrieval record lists",
                replace(
                    MINIMAL_CONFIGURATION,
                    "model = \"BSFG\"" => "model = \"BSFG\"\nmean_kinetic_energy_file = \"TKE_vs_A/2_B.Author_2001.dat\"",
                ),
            ),
            (
                "heavy mass range of one mass",
                replace(
                    MINIMAL_CONFIGURATION, "heavy_mass_max = 140" => "heavy_mass_max = 126"
                ),
            ),
            (
                "absent multiplicity subdirectory",
                replace(
                    MINIMAL_CONFIGURATION,
                    "subdirectory = \"datasets\"" => "subdirectory = \"absent\"",
                ),
            ),
            (
                "segment count out of range",
                replace(MINIMAL_CONFIGURATION, "max_segments = 3" => "max_segments = 99"),
            ),
            (
                "significant digits out of range",
                replace(
                    MINIMAL_CONFIGURATION, "significant_digits = 6" => "significant_digits = 0"
                ),
            ),
            (
                "window outside range",
                replace(
                    MINIMAL_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nrequired_windows = [[10, 20]]",
                ),
            ),
            (
                "malformed window",
                replace(
                    MINIMAL_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nrequired_windows = [[130]]",
                ),
            ),
            (
                "segment span out of range",
                replace(
                    MINIMAL_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nmin_segment_span = 0",
                ),
            ),
            (
                "coverage floor out of range",
                replace(
                    MINIMAL_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nmin_pair_coverage = 1.5",
                ),
            ),
            (
                "retrieval floor not a version",
                MINIMAL_CONFIGURATION * "\n[retrieval]\nmin_package_version = \"latest\"\n",
            ),
            (
                "unknown retrieval key",
                MINIMAL_CONFIGURATION * "\n[retrieval]\nversion = \"0.2.3\"\n",
            ),
            # A key the loader does not read is refused, not ignored: a retired or misspelt key
            # would otherwise leave the run silently on the default.
            (
                "retired key",
                replace(
                    MINIMAL_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nparsimony = 1.0",
                ),
            ),
            (
                "unknown key",
                replace(
                    MINIMAL_CONFIGURATION,
                    "significant_digits = 6" => "significant_digits = 6\nsubdirectory = \"x\"",
                ),
            ),
            ("unknown section", MINIMAL_CONFIGURATION * "\n[extra]\nkey = 1\n"),
        ]
        for (name, body) in cases
            @testset "$(name)" begin
                with_configuration(body) do path, directory
                    @test_throws ArgumentError load_configuration(
                        path; data_directory = directory
                    )
                end
            end
        end
    end

    @testset "pinning requires an even nucleus split at symmetry" begin
        body = replace(MINIMAL_CONFIGURATION, "target_A = 252" => "target_A = 251")
        with_configuration(body) do path, directory
            @test_throws ArgumentError load_configuration(path; data_directory = directory)
        end
        body = replace(
            replace(MINIMAL_CONFIGURATION, "target_A = 252" => "target_A = 251"),
            "max_segments = 3" => "max_segments = 3\npin_symmetric_split = false",
        )
        with_configuration(body) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test !has_symmetric_split(configuration)
        end

        # The fit is pinned at the first abscissa of the range, so a range starting above the
        # symmetric split would fix the ratio to one half where nothing says it is.
        body = replace(
            MINIMAL_CONFIGURATION,
            "heavy_mass_max = 140" => "heavy_mass_min = 130\nheavy_mass_max = 140",
        )
        with_configuration(body) do path, directory
            @test_throws ArgumentError load_configuration(path; data_directory = directory)
        end
        with_configuration(
            replace(body, "max_segments = 3" => "max_segments = 3\npin_symmetric_split = false")
        ) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test A_H_range(configuration) == 130:140
        end
    end

    @testset "the fissioning nucleus and the label are derived" begin
        # Neutron-induced fission adds one mass unit to the target; the label follows from the
        # target and the channel, so it cannot contradict the nuclide it names.
        body = replace(
            MINIMAL_CONFIGURATION,
            "target_A = 252\ntarget_Z = 98\nchannel = \"sf\"" => "target_A = 235\ntarget_Z = 92\nchannel = \"nth\"\nincident_energy = 2.53e-8",
        )
        body = replace(body, "heavy_mass_max = 140" => "heavy_mass_max = 160")
        with_configuration(body) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test configuration.system.compound == Nuclide(92, 236)
            @test reaction(configuration.system) == "n,f"
            @test system_label(configuration.system) == "U235_nth"
            @test configuration.system.incident_energy == 2.53e-8

            # Two incident energies of one target and channel are two systems, and must not
            # collide in the run identifier even though they share a label.
            other = joinpath(directory, "other.toml")
            write(other, replace(body, "incident_energy = 2.53e-8" => "incident_energy = 1.0"))
            @test run_identifier(configuration) !=
                run_identifier(load_configuration(other; data_directory = directory))
        end
    end

    @testset "a ⟨TKE⟩(A) dataset weights the charge-resolved inversion" begin
        body = replace(
            MINIMAL_CONFIGURATION,
            "model = \"BSFG\"" => "model = \"BSFG\"\nmean_kinetic_energy_file = \"TKE_vs_A/1_A.Author_2000.dat\"",
        )
        with_configuration(body) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test configuration.level_density.mean_kinetic_energy_file ==
                joinpath(directory, "TKE_vs_A", "1_A.Author_2000.dat")
            @test run_parameters(configuration)["TKE"] != "none"
            @test run_identifier(configuration) != run_identifier(
                load_configuration(
                    joinpath(write_beside(path, MINIMAL_CONFIGURATION));
                    data_directory = directory,
                ),
            )
        end
    end

    @testset "every result-changing key enters the run identifier" begin
        # The identifier is built from configuration keys through one exported table. A key
        # without an entry is an error rather than a silently invented token, and every entry is
        # used, so the table cannot carry a token for a key the identifier leaves out.
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            tokens = run_parameters(configuration)
            # The Gilbert-Cameron keys change a Gilbert-Cameron result only, the yield source
            # is a directory or one file, the yield coverage floor needs a yield, and the number
            # of correlogram lags is a token only away from its default.
            gilbert_cameron = Set(["sc", "def"])
            @test Set(keys(tokens)) == setdiff(
                Set(values(RUN_IDENTIFIER_ABBREVIATIONS)),
                gilbert_cameron,
                Set(["Yf", "Ycov", "lags"]),
            )
            other = write_beside(
                path, replace(MINIMAL_CONFIGURATION, "model = \"BSFG\"" => "model = \"GC\"")
            )
            gc_tokens = run_parameters(load_configuration(other; data_directory = directory))
            @test Set(keys(gc_tokens)) == setdiff(
                Set(values(RUN_IDENTIFIER_ABBREVIATIONS)), Set(["Yf", "Ycov", "lags"])
            )
            @test gc_tokens["sc"] == "gc1965"
            @test allunique(values(RUN_IDENTIFIER_ABBREVIATIONS))
            identifier = run_identifier(configuration)
            for token in keys(tokens)
                @test occursin("$(token)=", identifier)
            end
            # The system names the directory the identifier sits in, and is not a token.
            @test !occursin("Cf252", identifier)

            # The number of correlogram lags changes the uncertainty of the trend: a token, and
            # another run, where it is not the default.
            other = write_beside(
                path,
                replace(
                    MINIMAL_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nautocorrelation_lags = 2",
                ),
            )
            lagged = load_configuration(other; data_directory = directory)
            @test lagged.segments.autocorrelation_lags == 2
            @test configuration.segments.autocorrelation_lags == 4
            @test run_parameters(lagged)["lags"] == 2
            @test run_identifier(lagged) == "$(identifier)_lags=2" ||
                occursin("_lags=2_", run_identifier(lagged))

            # A list enters as a content hash, so two window sets are two runs.
            other = joinpath(directory, "windows.toml")
            write(
                other,
                replace(
                    MINIMAL_CONFIGURATION,
                    "max_segments = 3" => "max_segments = 3\nrequired_windows = [[130, 134]]",
                ),
            )
            @test run_identifier(load_configuration(other; data_directory = directory)) !=
                identifier

            # A path enters relative to the data directory, so the same inputs staged elsewhere
            # give the same identifier: a run made on another machine names the same inputs.
            elsewhere = joinpath(directory, "elsewhere")
            mkpath(joinpath(elsewhere, "datasets"))
            cp(joinpath(directory, "mass_excess.dat"), joinpath(elsewhere, "mass_excess.dat"))
            copy = joinpath(elsewhere, "configuration.toml")
            write(copy, MINIMAL_CONFIGURATION)
            @test run_identifier(load_configuration(copy; data_directory = elsewhere)) ==
                identifier
        end
    end

    @testset "the primary yield, alone or as the reference of a directory" begin
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            mkpath(joinpath(directory, "Y_vs_A"))
            write(joinpath(directory, "Y_vs_A", "1_A.Author_2000.dat"), "A Y\n130 5.0\n")
            # A tabulation that carries no accession, which an exclusion names by its label.
            write(joinpath(directory, "Y_vs_A", "2_B.Other_2001.dat"), "A Y\n130 4.0\n")
            file = write_beside(
                path,
                MINIMAL_CONFIGURATION *
                "\n[yield]\nmass_yield_file = \"Y_vs_A/1_A.Author_2000.dat\"\n",
            )
            configuration = load_configuration(file; data_directory = directory)
            @test configuration.yield_file ==
                joinpath(directory, "Y_vs_A", "1_A.Author_2000.dat")
            @test configuration.yield_directory === nothing
            tokens = run_parameters(configuration)
            @test haskey(tokens, "Yf") && !haskey(tokens, "Y")
            # An exclusion has nothing to act on where no directory is read, and is refused.
            excluding = write_beside(
                path,
                MINIMAL_CONFIGURATION *
                "\n[yield]\nmass_yield_file = \"Y_vs_A/1_A.Author_2000.dat\"\nexclude = [{ dataset = \"B. Other 2001\", reason = \"partial\" }]\n",
            )
            @test_throws "yield.exclude needs yield.subdirectory" load_configuration(
                excluding; data_directory = directory
            )
            # A malformed entry is reported for its form first.
            malformed = write_beside(
                path,
                MINIMAL_CONFIGURATION *
                "\n[yield]\nmass_yield_file = \"Y_vs_A/1_A.Author_2000.dat\"\nexclude = [{ dataset = \"B. Other 2001\" }]\n",
            )
            @test_throws "must be a table with a string `reason`" load_configuration(
                malformed; data_directory = directory
            )
            # A directory is averaged over with the primary as its coverage reference.
            both = write_beside(
                path,
                MINIMAL_CONFIGURATION *
                "\n[yield]\nsubdirectory = \"Y_vs_A\"\nmass_yield_file = \"Y_vs_A/1_A.Author_2000.dat\"\nexclude = [{ dataset = \"B. Other 2001\", reason = \"partial\" }]\n",
            )
            both = load_configuration(both; data_directory = directory)
            @test both.yield_directory == joinpath(directory, "Y_vs_A")
            @test both.yield_file == configuration.yield_file
            # One token for the directory, its reference and its exclusions.
            @test haskey(run_parameters(both), "Y") && !haskey(run_parameters(both), "Yf")
            unexcluded = write_beside(
                path,
                MINIMAL_CONFIGURATION *
                "\n[yield]\nsubdirectory = \"Y_vs_A\"\nmass_yield_file = \"Y_vs_A/1_A.Author_2000.dat\"\n",
            )
            @test run_parameters(both)["Y"] !=
                run_parameters(load_configuration(unexcluded; data_directory = directory))["Y"]
            # A directory without its reference has nothing to measure coverage against.
            directory_alone = write_beside(
                path, MINIMAL_CONFIGURATION * "\n[yield]\nsubdirectory = \"Y_vs_A\"\n"
            )
            @test_throws ArgumentError load_configuration(
                directory_alone; data_directory = directory
            )
            neither = write_beside(
                path, MINIMAL_CONFIGURATION * "\n[yield]\nsymmetrize = true\n"
            )
            @test_throws ArgumentError load_configuration(neither; data_directory = directory)
        end
    end

    @testset "two files of one accession are refused" begin
        # One measurement held twice, as a retrieval under another spelling of the author leaves
        # it: it would be read and pooled twice, and an exclusion could name only one.
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            for name in ("41739002_A.Set_1999.dat", "41739002_A.Sett_1999.dat")
                write(
                    joinpath(directory, "datasets", name), "A nu nu_uncertainty\n126 2.0 0.1\n"
                )
            end
            @test_throws "two files of the accession 41739002" load_configuration(
                path; data_directory = directory
            )
            @test_throws "two files of the accession 41739002" read_multiplicity_directory(
                joinpath(directory, "datasets")
            )
        end
        # A directory of yield distributions alike: one distribution would be averaged over
        # twice.
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            mkpath(joinpath(directory, "Y_vs_A"))
            write(joinpath(directory, "Y_vs_A", "1_A.Author_2000.dat"), "A Y\n130 5.0\n")
            for name in ("41739004_A.Set_1999.dat", "41739004_A.Sett_1999.dat")
                write(joinpath(directory, "Y_vs_A", name), "A Y Y_uncertainty\n134 2.0 0.1\n")
            end
            @test_throws "two files of the accession 41739004" read_mass_yield_directory(
                joinpath(directory, "Y_vs_A")
            )
            with_yields = write_beside(
                path,
                MINIMAL_CONFIGURATION *
                "\n[yield]\nsubdirectory = \"Y_vs_A\"\nmass_yield_file = \"Y_vs_A/1_A.Author_2000.dat\"\nexclude = [{ accession = \"41739004\", reason = \"x\" }]\n",
            )
            @test_throws "two files of the accession 41739004" load_configuration(
                with_yields; data_directory = directory
            )
        end
    end

    @testset "the retired coverage key names its replacements" begin
        body = replace(
            MINIMAL_CONFIGURATION,
            "max_segments = 3" => "max_segments = 3\nmin_dataset_coverage = 0.3",
        )
        with_configuration(body) do path, directory
            replacements = ["segments.min_pair_coverage", "yield.min_yield_coverage"]
            @test_throws replacements load_configuration(path; data_directory = directory)
        end
    end

    @testset "the yield coverage floor is its own key" begin
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            mkpath(joinpath(directory, "yields"))
            write(
                joinpath(directory, "yields", "10000001_A.Both_2000.dat"),
                "A Y Y_uncertainty\n112 1.0 0.1\n126 2.0 0.1\n140 1.0 0.1\n",
            )
            yield_section = "\n[yield]\nmass_yield_file = \"yields/10000001_A.Both_2000.dat\"\n"
            file = write_beside(
                path, MINIMAL_CONFIGURATION * yield_section * "min_yield_coverage = 0.6\n"
            )
            configuration = load_configuration(file; data_directory = directory)
            @test configuration.min_yield_coverage == 0.6
            @test configuration.segments.min_pair_coverage == 0.3
            @test run_parameters(configuration)["Ycov"] == 0.6
            @test run_parameters(configuration)["cov"] == 0.3
            out_of_range = write_beside(
                path, MINIMAL_CONFIGURATION * yield_section * "min_yield_coverage = 1.5\n"
            )
            @test_throws "yield.min_yield_coverage" load_configuration(
                out_of_range; data_directory = directory
            )
            # Without a yield distribution the floor gates nothing and carries no token.
            configuration = load_configuration(path; data_directory = directory)
            @test !haskey(run_parameters(configuration), "Ycov")
        end
    end

    @testset "retrieval records are held to a version floor" begin
        body = replace(
            MINIMAL_CONFIGURATION,
            "model = \"BSFG\"" => "model = \"BSFG\"\nmean_kinetic_energy_file = \"TKE_vs_A/1_A.Author_2000.dat\"",
        )
        accepted = "[[accepted]]\nfile = \"1_A.Author_2000.dat\"\nqualifiers = []\n"
        with_configuration(body) do path, directory
            record = joinpath(directory, "TKE_vs_A", "retrieval.toml")
            # The fixture's record was written by 0.2.3, the default floor.
            load_configuration(path; data_directory = directory)
            versions = check_retrieval_versions([joinpath(directory, "TKE_vs_A")], v"0.2.3")
            @test length(versions) == 1
            @test only(values(versions)) == v"0.2.3"

            write(record, accepted * "\n[run]\npackage_version = \"0.2.2\"\n")
            @test_throws "below retrieval.min_package_version" load_configuration(
                path; data_directory = directory
            )
            # A lower floor, set explicitly, admits it.
            lowered = write_beside(
                path, body * "\n[retrieval]\nmin_package_version = \"0.2.0\"\n"
            )
            configuration = load_configuration(lowered; data_directory = directory)
            @test configuration.min_retrieval_version == v"0.2.0"

            # A record that states no version cannot be held to the floor.
            write(record, accepted)
            @test_throws "states no [run] package_version" load_configuration(
                path; data_directory = directory
            )
            write(record, accepted * "\n[run]\npackage_revision = \"abc\"\n")
            @test_throws "states no [run] package_version" load_configuration(
                path; data_directory = directory
            )

            # Input directories written by different versions are admitted, and said to be.
            write(record, accepted * "\n[run]\npackage_version = \"0.2.3\"\n")
            write(
                joinpath(directory, "datasets", "retrieval.toml"),
                "[run]\npackage_version = \"0.2.4\"\n",
            )
            @test_logs (:warn, r"different versions") match_mode = :any load_configuration(
                path; data_directory = directory
            )
        end
        # A directory without a record, a tabulation not produced by the retrieval, is not judged.
        @test isempty(check_retrieval_versions([mktempdir()], v"9.9.9"))
    end

    @testset "the file names that carry the identifier fit a file system" begin
        # The longest the identifier gets: Gilbert-Cameron, a resonance energy, a ⟨TKE⟩(A)
        # dataset, exclusions, windows, a coverage floor with two decimals.
        body = replace(
            MINIMAL_CONFIGURATION,
            "target_A = 252\ntarget_Z = 98\nchannel = \"sf\"" => "target_A = 235\ntarget_Z = 92\nchannel = \"nres\"\nincident_energy = 5.8013e-4",
        )
        body = replace(
            body, "heavy_mass_max = 140" => "heavy_mass_min = 118\nheavy_mass_max = 160"
        )
        body = replace(
            body,
            "model = \"BSFG\"" => "model = \"GC\"\nratio_averaging = \"charge_resolved\"\nmean_kinetic_energy_file = \"TKE_vs_A/1_A.Author_2000.dat\"",
        )
        body = replace(
            body,
            "[multiplicity]" => "[multiplicity]\nexclude = [{ accession = \"41739002\", reason = \"b\" }]",
        )
        body = replace(
            body,
            "max_segments = 3" => "max_segments = 12\nrequired_windows = [[128, 132]]\nmin_pair_coverage = 0.35\nmin_points_per_segment = 10\nmin_segment_span = 10",
        )
        with_configuration(body) do path, directory
            write(
                joinpath(directory, "datasets", "41739002_A.Set_1999.dat"),
                "A nu nu_uncertainty\n126 2.0 0.1\n",
            )
            identifier = run_identifier(load_configuration(path; data_directory = directory))
            @test ncodeunits("total_average_R_T_$(identifier).csv") <= 255
        end
        # And every shipped configuration, under either level density model.
        if DATA_AVAILABLE
            mktempdir() do directory
                for system in SHIPPED_SYSTEMS,
                    model in ("BSFG", "GC"),
                    yields in (Dict(), Dict("subdirectory" => "$(system)/Y_vs_A"))

                    configuration = variant_configuration(
                        system,
                        directory;
                        set = Dict(
                            "level_density" => Dict("model" => model), "yield" => yields
                        ),
                    )
                    identifier = run_identifier(configuration)
                    @test ncodeunits("total_average_R_T_$(identifier).csv") <= 255
                end
            end
        else
            # Reported, not passed over: the measured input is not shipped.
            @test_skip DATA_AVAILABLE
        end
    end

    @testset "exclusions name a dataset by accession, with a written reason" begin
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            @test isempty(
                load_configuration(path; data_directory = directory).excluded_datasets
            )
        end
        excluding(clause) =
            replace(MINIMAL_CONFIGURATION, "[multiplicity]" => "[multiplicity]\n$(clause)")
        write_dataset(directory, file) =
            write(joinpath(directory, "datasets", file), "A nu nu_uncertainty\n126 2.0 0.1\n")

        # Excluding a measurement is a judgement; an unexplained one cannot be told from a slip.
        for clause in (
            "exclude = [{ accession = \"41739002\" }]",
            "exclude = [{ accession = \"41739002\", reason = \"  \" }]",
            "exclude = [\"41739002\"]",
            "exclude = [{ accession = \"41739002\", reason = \"x\" }, " *
            "{ accession = \"41739002\", reason = \"y\" }]",
            "exclude = [{ accession = \"41739002\", dataset = \"A. Set 1999\", reason = \"x\" }]",
            "exclude = [{ reason = \"x\" }]",
            "exclude = [{ accession = \"4173\", reason = \"x\" }]",
            "exclude = [{ accession = \"41739002\", reason = \"x\", note = \"y\" }]",
        )
            with_configuration(excluding(clause)) do path, directory
                write_dataset(directory, "41739002_A.Set_1999.dat")
                @test_throws ArgumentError load_configuration(path; data_directory = directory)
            end
        end

        # Two subentries of one author and year: each label carries its accession, and only
        # the accession names one of them whatever its neighbours are.
        by_accession = "exclude = [{ accession = \"41739002\", reason = \"too sparse\" }]"
        with_configuration(excluding(by_accession)) do path, directory
            write_dataset(directory, "41739002_A.Set_1999.dat")
            write_dataset(directory, "41739003_A.Set_1999.dat")
            configuration = @test_logs min_level = Base.CoreLogging.Warn load_configuration(
                path; data_directory = directory
            )
            @test configuration.excluded_datasets == Dict("41739002" => "too sparse")
        end
        by_label = "exclude = [{ dataset = \"A. Set 1999\", reason = \"too sparse\" }]"
        with_configuration(excluding(by_label)) do path, directory
            write_dataset(directory, "41739002_A.Set_1999.dat")
            write_dataset(directory, "41739003_A.Set_1999.dat")
            @test_throws "which is not among those of" load_configuration(
                path; data_directory = directory
            )
        end

        # An exclusion of a dataset the directory does not hold is refused.
        absent = "exclude = [{ accession = \"99999999\", reason = \"x\" }]"
        with_configuration(excluding(absent)) do path, directory
            write_dataset(directory, "41739002_A.Set_1999.dat")
            @test_throws "which no dataset of" load_configuration(
                path; data_directory = directory
            )
        end

        # A dataset that has an accession is excluded by it; its label is refused, with the
        # accession to write.
        with_configuration(excluding(by_label)) do path, directory
            write_dataset(directory, "41739002_A.Set_1999.dat")
            @test_throws "accession = \"41739002\"" load_configuration(
                path; data_directory = directory
            )
        end

        # A tabulation that carries no accession is named by its label.
        own = "exclude = [{ dataset = \"own\", reason = \"test\" }]"
        with_configuration(excluding(own)) do path, directory
            write_dataset(directory, "own.dat")
            configuration = @test_logs min_level = Base.CoreLogging.Warn load_configuration(
                path; data_directory = directory
            )
            @test configuration.excluded_datasets == Dict("own" => "test")
        end
    end

    @testset "the shipped configurations are named for the system they describe" begin
        # A configuration file cannot drift from the system it runs.
        directory = joinpath(pkgdir(FissionTemperatureRatio), "config")
        files = filter(endswith(".toml"), readdir(directory))
        @test !isempty(files)
        for file in files
            system = system_of(TOML.parsefile(joinpath(directory, file)))
            @test file == "$(system_label(system)).toml"
        end
    end
end
