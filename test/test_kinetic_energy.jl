# The ⟨TKE⟩(A) that weights the charge-resolved inversion, and the retrieval records every input
# dataset is traced through.

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

@testset "mean total kinetic energy" begin
    system = neutron_induced_fission(Nuclide(92, 233), 2.53e-8, "nth")

    @testset "carried onto the heavy-mass range" begin
        mktempdir() do directory
            path = joinpath(directory, "40112007_V.M.Surin_1972.dat")
            # Heavy wing 120:2:124 with a gap at 121 and 123, light-wing 113 = complement of 121,
            # and a light-wing 110 whose heavy complement 124 is also tabulated.
            write(
                path,
                "A TKE TKE_uncertainty\n120 163.0 0.0\n122 169.0 0.5\n124 174.0 1.0\n" *
                "113 166.0 0.0\n110 150.0 0.0\n",
            )
            energies = read_mean_kinetic_energy(path, system, 117:126)
            @test energies.label == "V.M. Surin 1972"
            @test sort(collect(keys(energies.values))) == collect(117:126)
            # The pair's TKE, tabulated at the light partner A₀ - A_H = 113.
            @test energies.values[121] == 166.0
            # Both wings of one split tabulated: the pair takes their mean.
            @test energies.values[124] ≈ (174.0 + 150.0) / 2
            @test energies.averaged == [124]
            # 123 lies between measured masses: interpolated.
            @test energies.values[123] ≈ (169.0 + 162.0) / 2
            @test energies.interpolated == [123]
            # Beyond the measured span: the nearest measured value.
            @test energies.values[117] == 163.0
            @test energies.values[126] == energies.values[124]
            @test energies.extrapolated == [117, 118, 119, 125, 126]
            @test energies.measured == [120, 121, 122, 124]
        end
    end

    @testset "offset from the energy standard" begin
        mktempdir() do directory
            path = joinpath(directory, "1_A.Author_2000.dat")
            write(path, "A TKE\n130 176.0\n140 170.0\n")
            energies = read_mean_kinetic_energy(path, system, 130:140)
            # Weights 1 and 3 at the two measured masses; the others carry no yield.
            yields = MassYield([130, 140], [1.0, 3.0], [missing, missing], "y", "")
            offset = mean_kinetic_energy_offset(energies, yields, system)
            @test offset.mean ≈ (176.0 + 3 * 170.0) / 4
            standard = recommended_mean_total_kinetic_energy(system)
            @test offset.standard == Measurements.value(standard) == 170.1
            @test offset.offset ≈ offset.mean - 170.1
            @test offset.masses == [130, 140]
            # A system without a recorded standard still gets its mean.
            resonance = neutron_induced_fission(Nuclide(92, 233), 1.0, "nres")
            @test mean_kinetic_energy_offset(energies, yields, resonance).offset === nothing
            @test_throws ArgumentError mean_kinetic_energy_offset(
                energies, MassYield([200], [1.0], [missing], "far", ""), system
            )
        end
    end

    @testset "rejects what is not a pre-neutron mean energy" begin
        mktempdir() do directory
            bad = joinpath(directory, "bad.dat")
            write(bad, "A TKE\n130 0.0\n")
            @test_throws ArgumentError read_mean_kinetic_energy(bad, system, 117:159)
            repeated = joinpath(directory, "repeated.dat")
            write(repeated, "A TKE\n130 170.0\n130 171.0\n")
            @test_throws ArgumentError read_mean_kinetic_energy(repeated, system, 117:159)
            # The one tabulated mass lies beyond the range: no heavy mass of it is measured.
            far = joinpath(directory, "far.dat")
            write(far, "A TKE\n200 170.0\n")
            @test_throws ArgumentError read_mean_kinetic_energy(far, system, 117:159)
        end
    end

    @testset "the retrieval record names the dataset and its qualifiers" begin
        mktempdir() do directory
            write(joinpath(directory, "retrieval.toml"), RETRIEVAL_RECORD)
            surin = joinpath(directory, "40112007_V.M.Surin_1972.dat")
            record = retrieval_record(surin)
            @test record !== nothing
            @test isempty(record.qualifiers)
            @test record.run["package_revision"] == "87a9533"
            @test record.entry["identifier"] == "40112007"

            flagged = retrieval_record(joinpath(directory, "21981008_P.Geltenbort_1985.dat"))
            @test flagged.qualifiers == ["MXW", "DERIV"]
            @test "DERIV" in FLAGGED_QUALIFIERS
            @test "MXW" ∉ FLAGGED_QUALIFIERS

            @test retrieval_record(joinpath(directory, "absent.dat")) === nothing
            # Neither dataset was interpolated onto the integers: full weight where pooled.
            @test pooling_weight(record) == 1.0
            @test pooling_weight(nothing) == 1.0
            interpolated = RetrievalRecord(
                "r",
                String[],
                Dict{String,Any}(
                    "mass_treatment" => "interpolated",
                    "mass_values_non_integer" => 52,
                    "rows_written" => 69,
                ),
                Dict{String,Any}(),
            )
            @test pooling_weight(interpolated) ≈ 52 / 69
            by_file = retrieval_qualifiers(directory)
            @test by_file["21981008_P.Geltenbort_1985.dat"] == ["MXW", "DERIV"]
            @test isempty(by_file["40112007_V.M.Surin_1972.dat"])
        end
    end
end
