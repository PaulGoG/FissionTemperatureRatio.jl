# The configuration parser must enforce exactly the constraints its comments document, and fail
# naming the offending key, so that a run cannot start from a configuration it cannot honour.

const MINIMAL_CONFIGURATION = """
[system]
label = "test"
A0 = 252
Z0 = 98

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
            @test configuration.output.subdirectory == "test"
        end
    end

    @testset "missing file" begin
        @test_throws ArgumentError load_configuration(joinpath(@__DIR__, "absent.toml"))
    end

    @testset "each constraint is enforced" begin
        cases = [
            ("missing section", replace(MINIMAL_CONFIGURATION, "[output]\ndigits = 6\n" => "")),
            ("missing key", replace(MINIMAL_CONFIGURATION, "A0 = 252\n" => "")),
            ("wrong type", replace(MINIMAL_CONFIGURATION, "A0 = 252" => "A0 = \"252\"")),
            (
                "even charges per mass",
                replace(
                    MINIMAL_CONFIGURATION, "charges_per_mass = 5" => "charges_per_mass = 4"
                ),
            ),
            ("charge above mass", replace(MINIMAL_CONFIGURATION, "Z0 = 98" => "Z0 = 300")),
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
        body = replace(MINIMAL_CONFIGURATION, "A0 = 252" => "A0 = 251")
        with_configuration(body) do path, directory
            @test_throws ArgumentError load_configuration(path; data_directory = directory)
        end
        body = replace(
            replace(MINIMAL_CONFIGURATION, "A0 = 252" => "A0 = 251"),
            "max_segments = 3" => "max_segments = 3\npin_symmetric_split = false",
        )
        with_configuration(body) do path, directory
            configuration = load_configuration(path; data_directory = directory)
            @test !has_symmetric_split(configuration)
        end
    end
end
