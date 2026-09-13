@testset "fragment mass yields" begin
    @testset "reading" begin
        mktempdir() do directory
            path = joinpath(directory, "40112004_V.M.Surin_1972.dat")
            write(path, "A yield erryield\n132 0.061 0.002\n130 0.058 0.001\n")
            data = read_yield(path)
            # Ascending in mass number regardless of the order in the file.
            @test data.A == [130, 132]
            @test data.Y ≈ [0.058, 0.061]
            @test data.σY ≈ [0.001, 0.002]
            @test mass_yield(data, 132) == (0.061, 0.002)
            @test ismissing(mass_yield(data, 131))
            @test length(data) == 2

            # A distribution quoting no uncertainties is carried, not discarded.
            two_column = joinpath(directory, "two_column.dat")
            write(two_column, "A yield\n130 0.058\n132 0.061\n")
            @test read_yield(two_column).σY == [0.0, 0.0]

            # Directory reading takes the label from the file name and ignores the run record.
            write(joinpath(directory, "retrieval.toml"), "[query]\nordinate = \"yield\"\n")
            rm(two_column)
            sets = read_yield_directory(directory)
            @test length(sets) == 1
            @test only(sets).label == "V.M. Surin 1972"
        end
    end

    @testset "rejects what is not data" begin
        mktempdir() do directory
            negative = joinpath(directory, "negative.dat")
            write(negative, "A yield\n130 -0.1\n")
            @test_throws ArgumentError read_yield(negative)

            repeated = joinpath(directory, "repeated.dat")
            write(repeated, "A yield\n130 0.05\n130 0.06\n")
            @test_throws ArgumentError read_yield(repeated)

            @test_throws ArgumentError read_yield(joinpath(directory, "absent.dat"))
            @test_throws ArgumentError read_yield_directory(joinpath(directory, "absent"))
        end
    end

    @testset "total average" begin
        curve = RatioCurve([130, 132], [1.2, 1.0], [0.0, 0.0], "example")
        exact = YieldData([130, 132], [1.0, 3.0], [0.0, 0.0], "exact", "")

        # ⟨R_T⟩ = (1·1.2 + 3·1.0)/4. With neither input carrying uncertainties the result is exact.
        mean, σ = total_average(curve, exact)
        @test mean ≈ 1.05
        @test σ == 0.0

        # The normalization of the yields cancels.
        for scale in (0.5, 100.0)
            scaled = YieldData(exact.A, scale .* exact.Y, exact.σY, "scaled", "")
            @test first(total_average(curve, scaled)) ≈ mean
        end

        # Ratio uncertainties propagate with the normalized yields as coefficients.
        from_ratio = total_average(
            RatioCurve([130, 132], [1.2, 1.0], [0.1, 0.1], "example"), exact
        )
        @test from_ratio[1] ≈ 1.05
        @test from_ratio[2] ≈ sqrt((1 / 4)^2 * 0.01 + (3 / 4)^2 * 0.01)

        # A data set quoting no uncertainties still gets one, from the yield distribution alone.
        # This is why the published table carries an uncertainty for such sets.
        from_yield = total_average(
            curve, YieldData([130, 132], [1.0, 3.0], [0.5, 0.5], "uncertain", "")
        )
        @test from_yield[1] ≈ 1.05
        @test from_yield[2] ≈ sqrt(((1.2 - 1.05) / 4)^2 * 0.25 + ((1.0 - 1.05) / 4)^2 * 0.25)
        @test from_yield[2] > 0

        # A yield distribution contributes nothing where the ratio sits at its own average.
        flat = RatioCurve([130, 132], [1.05, 1.05], [0.0, 0.0], "flat")
        @test last(
            total_average(flat, YieldData([130, 132], [1.0, 3.0], [0.5, 0.5], "y", ""))
        ) ≈ 0.0

        # Only the mass numbers the two share enter.
        partial = YieldData([132], [1.0], [0.0], "partial", "")
        @test first(total_average(curve, partial)) ≈ 1.0

        @test_throws ArgumentError total_average(
            curve, YieldData([150], [1.0], [0.0], "disjoint", "")
        )
        @test_throws ArgumentError total_average(
            curve, YieldData([130, 132], [0.0, 0.0], [0.0, 0.0], "empty", "")
        )
    end

    @testset "the total average is not the range mean" begin
        # The range mean weights every mass number equally, so the far-asymmetric tail, where the
        # yield is negligible, counts as much as the peak. The two must not be confused.
        A_H = collect(126:174)
        value = [a ≤ 140 ? 1.2 : 0.6 for a in A_H]
        curve = RatioCurve(A_H, value, fill(0.01, length(A_H)), "example")
        peaked = YieldData(
            A_H, [a ≤ 140 ? 1.0 : 1.0e-4 for a in A_H], zeros(length(A_H)), "y", ""
        )
        @test first(total_average(curve, peaked)) > first(weighted_mean(curve))
        @test first(total_average(curve, peaked)) ≈ 1.2 atol = 0.01
    end
end
