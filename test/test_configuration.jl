# The configuration parser must enforce exactly the constraints its comments document, and fail
# naming the offending key, so that a run cannot start from a configuration it cannot honour.

const MINIMAL_CONFIGURATION = """
[system]
target_A = 252
target_Z = 98
reaction = "0,f"

[fragmentation]
charges_per_mass = 5
A_H_max = 140

[level_density]
prescription = "BSFG"
mass_excess_file = "masses.ANA"

[multiplicity]
directory = "sets"

[segments]
max_segments = 3

[output]
digits = 6
"""

function with_configuration(f, body::AbstractString)
    mktempdir() do directory
        write(joinpath(directory, "masses.ANA"), "1 1 H 7289.0 0.0\n1 0 n 8071.0 0.0\n")
        mkpath(joinpath(directory, "sets"))
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
            @test A_H_min(configuration) == 126
            @test A_H_range(configuration) == 126:140
            @test has_symmetric_split(configuration)
            @test configuration.level_density.ratio_averaging isa RatioOfMeans
            @test configuration.segments.pin_symmetric_split
            @test configuration.output.subdirectory == "Cf252_0f"
        end
    end

    @testset "missing file" begin
        @test_throws ArgumentError load_configuration(joinpath(@__DIR__, "absent.toml"))
    end

    @testset "each constraint is enforced" begin
        cases = [
            ("missing section", replace(MINIMAL_CONFIGURATION, "[output]\ndigits = 6\n" => "")),
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
                "unknown reaction",
                replace(MINIMAL_CONFIGURATION, "reaction = \"0,f\"" => "reaction = \"p,f\""),
            ),
            (
                "energy given for spontaneous fission",
                replace(
                    MINIMAL_CONFIGURATION,
                    "reaction = \"0,f\"" => "reaction = \"0,f\"\nincident_energy = 1.0",
                ),
            ),
            (
                "range below symmetry",
                replace(MINIMAL_CONFIGURATION, "A_H_max = 140" => "A_H_max = 100"),
            ),
            (
                "range above A0",
                replace(MINIMAL_CONFIGURATION, "A_H_max = 140" => "A_H_max = 300"),
            ),
            (
                "unknown prescription",
                replace(
                    MINIMAL_CONFIGURATION, "prescription = \"BSFG\"" => "prescription = \"XYZ\""
                ),
            ),
            (
                "unknown averaging",
                MINIMAL_CONFIGURATION * "\n[level_density.extra]\n" |>
                s -> replace(
                    MINIMAL_CONFIGURATION,
                    "prescription = \"BSFG\"" => "prescription = \"BSFG\"\nratio_averaging = \"other\"",
                ),
            ),
            (
                "absent mass excess file",
                replace(
                    MINIMAL_CONFIGURATION,
                    "mass_excess_file = \"masses.ANA\"" => "mass_excess_file = \"absent.ANA\"",
                ),
            ),
            (
                "absent multiplicity directory",
                replace(
                    MINIMAL_CONFIGURATION, "directory = \"sets\"" => "directory = \"absent\""
                ),
            ),
            (
                "segment count out of range",
                replace(MINIMAL_CONFIGURATION, "max_segments = 3" => "max_segments = 99"),
            ),
            (
                "digits out of range",
                replace(MINIMAL_CONFIGURATION, "digits = 6" => "digits = 0"),
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

    @testset "pinning requires an even fissioning nucleus" begin
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
    end

    @testset "the fissioning nucleus and the label are derived" begin
        # Neutron-induced fission adds one mass unit to the target; the label follows from the
        # target and the reaction, so it cannot contradict the nuclide it names.
        body = replace(
            MINIMAL_CONFIGURATION,
            "target_A = 252\ntarget_Z = 98\nreaction = \"0,f\"" => "target_A = 235\ntarget_Z = 92\nreaction = \"n,f\"\nincident_energy = 2.53e-8",
        )
        with_configuration(body) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test configuration.system.A₀ == 236
            @test configuration.system.Z₀ == 92
            @test configuration.system.label == "U235_nf"
            @test configuration.system.incident_energy == 2.53e-8

            # Two incident energies of one target and reaction are two systems, and must not
            # collide in the run identifier even though they share a label.
            other = joinpath(directory, "other.toml")
            write(other, replace(body, "incident_energy = 2.53e-8" => "incident_energy = 1.0"))
            @test run_identifier(configuration) !=
                run_identifier(load_configuration(other; data_directory = directory))
        end

        @test case_label(252, 98, "0,f") == "Cf252_0f"
        @test case_label(235, 92, "n,f") == "U235_nf"
        @test element_symbol(98) == "Cf"
        @test_throws ArgumentError element_symbol(0)
    end

    @testset "excluded data sets need a written reason" begin
        with_configuration(MINIMAL_CONFIGURATION) do path, directory
            @test isempty(load_configuration(path; data_directory = directory).excluded_sets)
        end

        # Excluding a measurement is a judgement; an unexplained one cannot be told from a slip.
        for clause in (
            "exclude = [{ set = \"A\" }]",
            "exclude = [{ set = \"A\", reason = \"  \" }]",
            "exclude = [\"A\"]",
            "exclude = [{ set = \"A\", reason = \"x\" }, { set = \"A\", reason = \"y\" }]",
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
            "[multiplicity]" => "[multiplicity]\nexclude = [{ set = \"A. Set 1999\", reason = \"too sparse\" }]",
        )
        with_configuration(body) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test configuration.excluded_sets["A. Set 1999"] == "too sparse"
        end
    end
end
