# The retrieval records every input dataset is traced through, read once and strictly.

const RETRIEVAL_RECORD = """
[[accepted]]
file = "40112007_V.M.Surin_1972.dat"
identifier = "40112007"
qualifiers = []
reaction_code = "92-U-233(N,F)MASS,PRE,KE,LF+HF"
[[accepted]]
file = "21981008_P.Geltenbort_1985.dat"
identifier = "21981008"
qualifiers = ["MXW: Maxwellian-averaged", "DERIV: derived from other data rather than measured"]
reaction_code = "92-U-233(N,F)MASS,PRE,KE,LF+HF,MXW,DERIV"
[[rejected]]
identifier = "416110032"
reason = "missing required tag \\"MASS\\""

[run]
package_revision = "87a9533"
package_version = "0.1.0"
"""

# A record of one accepted entry, the file `20000001_A.First_1979.dat`, with `keys` added to it.
function single_entry_record(keys::AbstractString)
    return """
        [[accepted]]
        file = "20000001_A.First_1979.dat"
        identifier = "20000001"
        $(keys)

        [run]
        package_version = "0.7.0"
        """
end

@testset "retrieval records" begin
    @testset "the entry names the dataset and its qualifiers" begin
        mktempdir() do directory
            write(joinpath(directory, "retrieval.toml"), RETRIEVAL_RECORD)
            records = retrieval_records(directory)
            @test sort!(collect(keys(records))) ==
                ["21981008_P.Geltenbort_1985.dat", "40112007_V.M.Surin_1972.dat"]
            surin = records["40112007_V.M.Surin_1972.dat"]
            @test surin.accession == "40112007"
            @test isempty(surin.qualifiers)
            @test surin.run["package_revision"] == "87a9533"
            @test surin.entry["reaction_code"] == "92-U-233(N,F)MASS,PRE,KE,LF+HF"
            # Neither dataset was interpolated onto the integers: full weight where pooled.
            @test surin.measured_fraction == 1.0
            @test surin.pair_sum === nothing
            @test isempty(surin.correlated_with) && isempty(surin.correlation_relation)
            @test records["21981008_P.Geltenbort_1985.dat"].qualifiers == ["MXW", "DERIV"]

            @test retrieval_record(joinpath(directory, "40112007_V.M.Surin_1972.dat")).path ==
                surin.path
            @test retrieval_record(joinpath(directory, "absent.dat")) === nothing
            @test isempty(retrieval_records(mktempdir()))
        end
    end

    @testset "measured fraction, experiment and pair-sum scale" begin
        mktempdir() do directory
            path = joinpath(directory, "retrieval.toml")
            write(path, single_entry_record("""
                        mass_treatment = "interpolated"
                        mass_values_non_integer = 52
                        rows_written = 69
                        correlated_with = ["20000002"]
                        correlation_relation = "alternative_analysis"
                        pair_sum_deviation = -0.1071
                        pair_sum_deviation_uncertainty = 0.0133
                        scale_consistent = false
                        """))
            record = only(values(retrieval_records(directory)))
            @test record.measured_fraction ≈ 52 / 69
            # More raw masses than rows: a measurement behind every row, and no more weight.
            write(path, single_entry_record("""
                        mass_treatment = "interpolated"
                        mass_values_non_integer = 12
                        rows_written = 8
                        """))
            @test only(values(retrieval_records(directory))).measured_fraction == 1.0
            write(path, single_entry_record("""
                        mass_treatment = "interpolated"
                        mass_values_non_integer = 52
                        rows_written = 69
                        correlated_with = ["20000002"]
                        correlation_relation = "alternative_analysis"
                        pair_sum_deviation = -0.1071
                        pair_sum_deviation_uncertainty = 0.0133
                        scale_consistent = false
                        """))
            record = only(values(retrieval_records(directory)))
            @test record.correlated_with == ["20000002"]
            @test record.correlation_relation == "alternative_analysis"
            @test record.pair_sum ==
                FissionTemperatureRatio.PairSumScale(-0.1071, 0.0133, false)

            # Without a stated uncertainty the record forms no verdict on the scale.
            write(path, single_entry_record("pair_sum_deviation = -0.0069"))
            unweighed = only(values(retrieval_records(directory))).pair_sum
            @test ismissing(unweighed.uncertainty) && ismissing(unweighed.consistent)
        end
    end

    @testset "a key of the wrong type is a defect of the record" begin
        mktempdir() do directory
            path = joinpath(directory, "retrieval.toml")
            for (keys, named) in (
                # A count that is not an integer would pool the dataset at full weight.
                (
                    "mass_treatment = \"interpolated\"\nmass_values_non_integer = 52.5\n" *
                    "rows_written = 69",
                    "mass_values_non_integer",
                ),
                (
                    "mass_treatment = \"interpolated\"\nmass_values_non_integer = true\n" *
                    "rows_written = 69",
                    "mass_values_non_integer",
                ),
                (
                    "mass_treatment = \"interpolated\"\nmass_values_non_integer = 0\n" *
                    "rows_written = 69",
                    "mass_values_non_integer",
                ),
                # A string where an array stands would pool one experiment twice.
                ("correlated_with = \"20000002\"", "correlated_with"),
                ("correlated_with = [\"two\"]", "correlated_with"),
                ("correlation_relation = 1", "correlation_relation"),
                ("qualifiers = \"MXW\"", "qualifiers"),
                ("pair_sum_deviation = \"-0.1\"", "pair_sum_deviation"),
                ("pair_sum_deviation = 0.1\nscale_consistent = 1", "scale_consistent"),
            )
                write(path, single_entry_record(keys))
                error = try
                    retrieval_records(directory)
                    nothing
                catch exception
                    exception
                end
                @test error isa ArgumentError
                @test occursin(named, sprint(showerror, error))
                @test occursin("20000001_A.First_1979.dat", sprint(showerror, error))
            end

            write(
                path,
                replace(
                    single_entry_record(""),
                    "identifier = \"20000001\"" => "identifier = \"2000\"",
                ),
            )
            @test_throws ArgumentError retrieval_records(directory)
            # One file listed by two records: which entry holds would depend on the order read.
            write(path, single_entry_record(""))
            write(joinpath(directory, "retrieval_second.toml"), single_entry_record(""))
            @test_throws ArgumentError retrieval_records(directory)
        end
    end

    @testset "the multiplicity and ⟨TKE⟩ readers carry the record" begin
        mktempdir() do directory
            write(joinpath(directory, "retrieval.toml"), single_entry_record(""))
            path = joinpath(directory, "20000001_A.First_1979.dat")
            write(path, "A nu nu_uncertainty\n126 2.0 0.1\n")
            @test read_multiplicity(path).record.accession == "20000001"
            @test only(read_multiplicity_directory(directory)).record.accession == "20000001"
            # A tabulation no record lists carries none.
            @test read_multiplicity(path; record = nothing).record === nothing
        end
    end
end
