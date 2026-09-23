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
mass_excess_file = "mass_excess.dat"

[multiplicity]
subdirectory = "datasets"

[segments]
max_segments = 3

[output]
significant_digits = 6
"""

function with_configuration(f, body::AbstractString)
    mktempdir() do directory
        write(joinpath(directory, "mass_excess.dat"), "1 1 H 7289.0 0.0\n1 0 n 8071.0 0.0\n")
        mkpath(joinpath(directory, "datasets"))
        path = joinpath(directory, "configuration.toml")
        write(path, body)
        return f(path, directory)
    end
end

@testset "configuration" begin
    @testset "a valid configuration parses" begin
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test configuration.system.A₀ == 252
            @test configuration.system.channel == "sf"
            @test configuration.system.reaction == "0,f"
            # The smallest heavy mass defaults to the symmetric split.
            @test configuration.fragmentation.heavy_mass_min == 126
            @test A_H_range(configuration) == 126:140
            @test has_symmetric_split(configuration)
            @test configuration.level_density.ratio_averaging isa RatioOfMeans
            @test configuration.fragmentation.fallback_charge_dispersion == 0.6
            @test configuration.segments.pin_symmetric_split
            @test configuration.segments.min_segment_span == 3
            @test configuration.segments.min_dataset_coverage == 0.3
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
                "Gilbert-Cameron without its shell corrections",
                replace(MINIMAL_CONFIGURATION, "model = \"BSFG\"" => "model = \"GC\""),
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
                    "mass_excess_file = \"mass_excess.dat\"" => "mass_excess_file = \"absent.dat\"",
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
                    "max_segments = 3" => "max_segments = 3\nmin_dataset_coverage = 1.5",
                ),
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
            @test configuration.system.A₀ == 236
            @test configuration.system.Z₀ == 92
            @test configuration.system.reaction == "n,f"
            @test configuration.system.label == "U235_nth"
            @test configuration.system.incident_energy == 2.53e-8

            # Two incident energies of one target and channel are two systems, and must not
            # collide in the run identifier even though they share a label.
            other = joinpath(directory, "other.toml")
            write(other, replace(body, "incident_energy = 2.53e-8" => "incident_energy = 1.0"))
            @test run_identifier(configuration) !=
                run_identifier(load_configuration(other; data_directory = directory))
        end

        @test system_label(252, 98, "sf") == "Cf252_sf"
        @test system_label(235, 92, "nth") == "U235_nth"
        @test system_label(235, 92, "nres") == "U235_nres"
        @test element_symbol(98) == "Cf"
        @test_throws ArgumentError element_symbol(0)
    end

    @testset "the typeset notation is a separate name from the token" begin
        # One name may not mean both the path token and the figure label.
        @test system_notation(252, 98, "sf") == "²⁵²Cf(sf)"
        @test system_notation(233, 92, "nth") == "²³³U(nth,f)"
        @test system_notation(235, 92, "nres") == "²³⁵U(nres,f)"
        @test system_notation(252, 98, "sf") != system_label(252, 98, "sf")
        @test_throws ArgumentError system_notation(252, 98, "0,f")
    end

    @testset "every result-changing key enters the run identifier" begin
        # The identifier is built from configuration keys through one exported table. A key
        # without an entry is an error rather than a silently invented token, and every entry is
        # used, so the table cannot carry a token for a key the identifier leaves out.
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            tokens = run_parameters(configuration)
            @test Set(keys(tokens)) == Set(values(RUN_IDENTIFIER_ABBREVIATIONS))
            @test allunique(values(RUN_IDENTIFIER_ABBREVIATIONS))
            identifier = run_identifier(configuration)
            for token in keys(tokens)
                @test occursin("$(token)=", identifier)
            end
            # The system names the directory the identifier sits in, and is not a token.
            @test !occursin("Cf252", identifier)

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

    @testset "excluded datasets need a written reason" begin
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            @test isempty(
                load_configuration(path; data_directory = directory).excluded_datasets
            )
        end

        # Excluding a measurement is a judgement; an unexplained one cannot be told from a slip.
        for clause in (
            "exclude = [{ dataset = \"A\" }]",
            "exclude = [{ dataset = \"A\", reason = \"  \" }]",
            "exclude = [\"A\"]",
            "exclude = [{ dataset = \"A\", reason = \"x\" }, { dataset = \"A\", reason = \"y\" }]",
        )
            body = replace(
                MINIMAL_CONFIGURATION, "[multiplicity]" => "[multiplicity]\n$(clause)"
            )
            with_configuration(body) do path, directory
                @test_throws ArgumentError load_configuration(path; data_directory = directory)
            end
        end

        body = replace(
            MINIMAL_CONFIGURATION,
            "[multiplicity]" => "[multiplicity]\nexclude = [{ dataset = \"A. Set 1999\", reason = \"too sparse\" }]",
        )
        with_configuration(body) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test configuration.excluded_datasets["A. Set 1999"] == "too sparse"
        end
    end

    @testset "the shipped configurations are named for the system they describe" begin
        # A configuration file cannot drift from the system it runs.
        directory = joinpath(pkgdir(FissionTemperatureRatio), "config")
        files = filter(endswith(".toml"), readdir(directory))
        @test !isempty(files)
        for file in files
            document = TOML.parsefile(joinpath(directory, file))
            label = system_label(
                document["system"]["target_A"],
                document["system"]["target_Z"],
                document["system"]["channel"],
            )
            @test file == "$(label).toml"
        end
    end
end
