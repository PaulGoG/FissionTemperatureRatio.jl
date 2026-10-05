
@testset "multiplicity ratio" begin
    A₀ = 252
    data = Multiplicity([120, 126, 132], [1.0, 2.0, 3.0], [0.1, 0.2, 0.3], "synthetic", "-")

    @testset "definition and symmetry" begin
        ratio = multiplicity_ratio(data, A₀, 126:132)
        # The symmetric split pairs a fragment with itself, so the ratio is one half exactly.
        index = findfirst(==(126), ratio.A_H)
        @test ratio.ratio[index] ≈ 0.5
        # ν(132) = 3 against its complement ν(120) = 1.
        index = findfirst(==(132), ratio.A_H)
        @test ratio.ratio[index] ≈ 3 / 4
    end

    @testset "uncertainty propagation" begin
        ratio = multiplicity_ratio(data, A₀, 132:132)
        ν_H, σ_H, ν_L, σ_L = 3.0, 0.3, 1.0, 0.1
        expected = sqrt((ν_L * σ_H)^2 + (ν_H * σ_L)^2) / (ν_L + ν_H)^2
        @test only(ratio.σ) ≈ expected
    end

    @testset "unpaired and unphysical points are omitted, never clipped" begin
        unpaired = Multiplicity([132], [3.0], [missing], "unpaired", "-")
        @test isempty(multiplicity_ratio(unpaired, A₀, 126:140))
        # A vanishing light multiplicity puts the ratio at the singular endpoint of the
        # temperature ratio relation; the point must be dropped rather than pushed to 1.
        degenerate = Multiplicity([120, 132], [0.0, 3.0], [missing, missing], "degenerate", "-")
        ratio = multiplicity_ratio(degenerate, A₀, 132:132)
        @test isempty(ratio)
    end

    @testset "reader stores an absent uncertainty as missing, keeping the point" begin
        mktempdir() do directory
            path = joinpath(directory, "set.dat")
            write(path, "A nu nu_uncertainty\n126 2.0 0.1\n132 3.0\n")
            data = read_multiplicity(path; label = "set")
            @test length(data) == 2
            @test ismissing(data.σν[2])
            # A retrieval that writes a blank uncertainty as 0.0 has quoted none: no measured
            # multiplicity is exact.
            zero = joinpath(directory, "zero.dat")
            write(zero, "A nu nu_uncertainty\n126 2.0 0.1\n132 3.0 0.0\n")
            @test ismissing(read_multiplicity(zero; label = "zero").σν[2])
        end
    end

    @testset "reader refuses a field it cannot read, naming the file and the line" begin
        mktempdir() do directory
            file(name, body) = (path = joinpath(directory, name); write(path, body); path)
            header = "A nu nu_uncertainty\n"
            # A multiplicity that is not a number: an ArgumentError, as the reader promises.
            value = file("value.dat", header * "100 1.5 0.1\n101 abc 0.1\n102 1.7 0.1\n")
            @test_throws ArgumentError read_multiplicity(value)
            @test_throws "value.dat, line 3" read_multiplicity(value)
            @test_throws "the multiplicity \"abc\" is not numeric" read_multiplicity(value)
            # One uncertainty that is not a number must not pass for a dataset that quotes none.
            uncertainty = file("uncertainty.dat", header * "100 1.5 0.1\n101 1.6 n/a\n")
            @test_throws ArgumentError read_multiplicity(uncertainty)
            @test_throws "uncertainty.dat, line 3" read_multiplicity(uncertainty)
            @test_throws "the uncertainty \"n/a\" is not numeric" read_multiplicity(uncertainty)
            # A mass number that is not a number, and a line that stops after it.
            mass = file("mass.dat", header * "x 1.5 0.1\n")
            @test_throws "mass.dat, line 2" read_multiplicity(mass)
            short = file("short.dat", header * "100 1.5 0.1\n101\n")
            @test_throws "short.dat, line 3: fewer than two fields" read_multiplicity(short)
            negative = file("negative.dat", header * "100 -1.5 0.1\n")
            @test_throws "negative.dat, line 2" read_multiplicity(negative)
            infinite = file("infinite.dat", header * "100 Inf 0.1\n")
            @test_throws ArgumentError read_multiplicity(infinite)
            empty = file("empty.dat", header)
            @test_throws "contains no usable rows" read_multiplicity(empty)
            twice = file("twice.dat", header * "100 1.5 0.1\n100 1.6 0.1\n")
            @test_throws "repeated mass numbers" read_multiplicity(twice)

            # What is read: the quoted uncertainties stay with their points; zero, a negative
            # number and NaN are "none quoted"; blank lines, tabs and further fields are passed.
            mixed = file(
                "mixed.dat",
                header *
                "102 1.7 0.2 extra\n\n100\t1.5\t0.1\n101 1.6 NaN\n103 1.8 -1\n104 1.9 0\n105 2.0\n",
            )
            data = read_multiplicity(mixed; label = "mixed")
            @test data.A == [100, 101, 102, 103, 104, 105]
            @test data.ν == [1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
            @test isequal(data.σν, [0.1, missing, 0.2, missing, missing, missing])
            @test data.label == "mixed"
        end
    end

    @testset "pooled datasets are weighted by their measured points" begin
        a = RatioCurve([130, 131], [0.30, 0.40], [0.01, 0.01], "a")
        b = RatioCurve([130], [0.31], [0.01], "b")
        # Values consistent within their uncertainties, so no between-dataset variance enters and
        # the combination is the weighted mean.
        equal = consensus([a, b])
        @test equal.ratio[1] ≈ (0.30 + 0.31) / 2
        halved = consensus([a, b]; weights = [1.0, 0.5])
        @test halved.ratio[1] ≈ (0.30 + 0.5 * 0.31) / 1.5
        @test halved.σ[1] ≈ 0.01 / sqrt(1.5)
        # Alone at a mass number, a curve passes through with its uncertainty inflated by 1/√f.
        alone = consensus([a, b]; weights = [0.25, 1.0])
        @test alone.ratio[2] == 0.40
        @test alone.σ[2] ≈ 0.01 / sqrt(0.25)
        @test_throws DimensionMismatch consensus([a, b]; weights = [1.0])
        @test_throws ArgumentError consensus([a, b]; weights = [1.0, 0.0])
    end
end
