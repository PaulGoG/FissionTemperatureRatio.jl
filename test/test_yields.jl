@testset "fragment mass yields" begin
    # The reader of one distribution is FissionFragmentsDomain's; the directory is read here.
    @testset "reading a directory" begin
        mktempdir() do directory
            write(
                joinpath(directory, "40112004_V.M.Surin_1972.dat"),
                "A Y Y_uncertainty\n132 0.061 0.002\n130 0.058 0.001\n",
            )
            # The run record beside the data is not data.
            write(joinpath(directory, "retrieval.toml"), "[query]\nordinate = \"yield\"\n")
            sets = read_mass_yield_directory(directory)
            @test length(sets) == 1
            data = only(sets)
            # The label from the file name, as for the multiplicity datasets.
            @test data.label == "V.M. Surin 1972"
            @test data.A == [130, 132]
            @test mass_yield(data, 132) == (0.061, 0.002)
            @test mass_yield(data, 131) === nothing

            # Two distributions of one author and year keep their identifiers: total averages are
            # keyed by label, and a shared one would let the second replace the first.
            for identifier in ("23717003", "23717005")
                write(
                    joinpath(directory, "$(identifier)_G.Barreau_1985.dat"),
                    "A Y Y_uncertainty\n140 1.0 0.1\n",
                )
            end
            @test [y.label for y in read_mass_yield_directory(directory)] == [
                "G. Barreau 1985 (23717003)", "G. Barreau 1985 (23717005)", "V.M. Surin 1972"
            ]
            for identifier in ("23717003", "23717005")
                rm(joinpath(directory, "$(identifier)_G.Barreau_1985.dat"))
            end

            @test_throws ArgumentError read_mass_yield_directory(joinpath(directory, "absent"))
            empty = joinpath(directory, "empty")
            mkpath(empty)
            @test_throws ArgumentError read_mass_yield_directory(empty)
        end
    end

    @testset "symmetrized to the pre-neutron identity" begin
        # A₀ = 236: 118 is its own complement, 130 ↔ 106 both measured, 140 has no partner and
        # stands for its complement 96.
        yields = MassYield(
            [106, 118, 130, 140],
            [2.0, 0.5, 4.0, 6.0],
            [0.3, missing, 0.4, 0.2],
            "y",
            "source.dat",
        )
        symmetric = symmetrized_mass_yield(yields, 236)
        @test symmetric.A == [96, 106, 118, 130, 140]
        @test symmetric.Y ≈ [6.0, 3.0, 0.5, 3.0, 6.0]
        @test symmetric.σY[2] ≈ sqrt(0.3^2 + 0.4^2) / 2
        @test symmetric.σY[4] ≈ symmetric.σY[2]
        @test ismissing(symmetric.σY[3])
        @test symmetric.σY[1] == symmetric.σY[5] == 0.2
        # The identity holds at every mass, and imposing it twice changes nothing.
        @test all(
            isequal(mass_yield(symmetric, A), mass_yield(symmetric, 236 - A)) for
            A in symmetric.A
        )
        twice = symmetrized_mass_yield(symmetric, 236)
        @test twice.A == symmetric.A && twice.Y ≈ symmetric.Y
        @test (symmetric.label, symmetric.source) == (yields.label, yields.source)
        # A distribution measured on the light wing alone gives the heavy-fragment yields.
        light = MassYield([100, 106, 110], [1.0, 2.0, 1.5], [0.1, 0.2, missing], "light", "")
        heavy = symmetrized_mass_yield(light, 236)
        @test isequal(
            [mass_yield(heavy, A) for A in (126, 130, 136)],
            [(1.5, missing), (2.0, 0.2), (1.0, 0.1)],
        )
        # One quoted uncertainty of a pair: the other contributes nothing.
        partial = MassYield([106, 130], [2.0, 4.0], [missing, 0.4], "p", "")
        @test symmetrized_mass_yield(partial, 236).σY == [0.2, 0.2]
    end

    @testset "total average" begin
        curve = RatioCurve([130, 132], [1.2, 1.0], [missing, missing], "example")
        exact = MassYield([130, 132], [1.0, 3.0], [missing, missing], "exact", "")

        # ⟨R_T⟩ = (1·1.2 + 3·1.0)/4. With neither input carrying uncertainties the result is exact.
        mean, σ = total_average(curve, exact)
        @test mean ≈ 1.05
        @test σ == 0.0

        # The normalization of the yields cancels.
        for scale in (0.5, 100.0)
            scaled = MassYield(exact.A, scale .* exact.Y, exact.σY, "scaled", "")
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
            curve, MassYield([130, 132], [1.0, 3.0], [0.5, 0.5], "uncertain", "")
        )
        @test from_yield[1] ≈ 1.05
        @test from_yield[2] ≈ sqrt(((1.2 - 1.05) / 4)^2 * 0.25 + ((1.0 - 1.05) / 4)^2 * 0.25)
        @test from_yield[2] > 0

        # A yield distribution contributes nothing where the ratio sits at its own average.
        flat = RatioCurve([130, 132], [1.05, 1.05], [missing, missing], "flat")
        @test last(
            total_average(flat, MassYield([130, 132], [1.0, 3.0], [0.5, 0.5], "y", ""))
        ) ≈ 0.0

        # Only the mass numbers the two share enter.
        partial = MassYield([132], [1.0], [missing], "partial", "")
        @test first(total_average(curve, partial)) ≈ 1.0

        @test_throws ArgumentError total_average(
            curve, MassYield([150], [1.0], [missing], "disjoint", "")
        )
        @test_throws ArgumentError total_average(
            curve, MassYield([130, 132], [0.0, 0.0], [missing, missing], "empty", "")
        )
    end

    @testset "coverage of the fragmentation range and fraction of the yield" begin
        # Coverage is measured in the reference's yield, not the set's own: three tail masses
        # normalised to one another cover what the reference holds there.
        reference = MassYield(
            collect(126:174),
            [A < 160 ? 2.0 : 0.5 for A in 126:174],
            fill(missing, 49),
            "ref",
            "",
        )
        tail = MassYield([170, 171, 172], [1.0, 1.0, 1.0], fill(missing, 3), "tail", "")
        @test mass_yield_coverage(tail, reference, 126:174) ≈ 1.5 / (34 * 2.0 + 15 * 0.5)
        @test mass_yield_coverage(tail, reference, 170:172) == 1.0
        @test mass_yield_coverage(reference, reference, 126:174) == 1.0
        @test_throws ArgumentError mass_yield_coverage(tail, reference, 180:190)

        yields = MassYield([130, 132, 134], [1.0, 2.0, 1.0], fill(missing, 3), "y", "")
        @test yield_fraction(
            RatioCurve([130, 132], [1.0, 1.0], [missing, missing], "c"), yields, 126:140
        ) ≈ 0.75
        @test yield_fraction(
            RatioCurve(collect(126:140), ones(15), fill(missing, 15), "c"), yields, 126:140
        ) == 1.0
        # Yield outside the range is not counted against the curve.
        @test yield_fraction(
            RatioCurve([130, 132], [1.0, 1.0], [missing, missing], "c"), yields, 130:132
        ) == 1.0
        @test_throws ArgumentError yield_fraction(
            RatioCurve([130], [1.0], [missing], "c"), yields, 140:150
        )
    end

    @testset "the total average is not the range mean" begin
        # The range mean weights every mass number equally, so the far-asymmetric tail, where the
        # yield is negligible, counts as much as the peak. The two must not be confused.
        A_H = collect(126:174)
        value = [a ≤ 140 ? 1.2 : 0.6 for a in A_H]
        curve = RatioCurve(A_H, value, fill(0.01, length(A_H)), "example")
        peaked = MassYield(
            A_H, [a ≤ 140 ? 1.0 : 1.0e-4 for a in A_H], zeros(length(A_H)), "y", ""
        )
        @test first(total_average(curve, peaked)) > sum(curve.ratio) / length(curve)
        @test first(total_average(curve, peaked)) ≈ 1.2 atol = 0.01
    end
end
