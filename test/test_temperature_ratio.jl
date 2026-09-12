
@testset "temperature ratio" begin
    @testset "equal sharing gives equal temperatures" begin
        # r_ν = 1/2 with R_a = 1 means the pair shares its excitation energy equally between two
        # fragments of equal level density parameter, so their temperatures coincide.
        curve = RatioCurve([126], [0.5], [0.0], "synthetic")
        R_T = temperature_ratio(curve, Dict(126 => 1.0))
        @test only(R_T.value) ≈ 1.0
    end

    @testset "relation is the inverse of the excitation energy partition" begin
        R_a = 1.3
        for r in (0.2, 0.35, 0.5, 0.7)
            curve = RatioCurve([130], [r], [0.0], "synthetic")
            R_T = only(temperature_ratio(curve, Dict(130 => R_a)).value)
            # E*_H/TXE = 1/(1 + R_a R_T²) must return the ratio it was built from.
            @test 1 / (1 + R_a * R_T^2) ≈ r
        end
    end

    @testset "uncertainty propagation" begin
        r, σ_r, R_a = 0.3, 0.02, 1.2
        curve = RatioCurve([130], [r], [σ_r], "synthetic")
        result = temperature_ratio(curve, Dict(130 => R_a))
        R_T = only(result.value)
        @test only(result.σ) ≈ σ_r / (2 * R_T * R_a * r^2)
    end

    @testset "singular and absent points are omitted" begin
        curve = RatioCurve([126, 130], [0.5, 0.4], [0.0, 0.0], "synthetic")
        @test length(temperature_ratio(curve, Dict(126 => 1.0))) == 1
    end

    @testset "averaging over the charge distribution" begin
        A₀, Z₀ = 252, 98
        domain = fragmentation_domain(A₀, Z₀, 126:132, 5, FLAT_CHARGES)
        prescription = BackShiftedFermiGas()
        by_means = level_density_ratio(
            RatioOfMeans(), prescription, A₀, Z₀, domain, TEST_MASSES
        )
        by_ratios = level_density_ratio(
            MeanOfRatios(), prescription, A₀, Z₀, domain, TEST_MASSES
        )

        # At the symmetric split the two fragments are the same nuclide, so the exact ratio is one.
        # Averaging the parameters and then dividing reproduces that identity exactly, because the
        # retained charges are invariant under Z -> Z₀ - Z and the two averages run over the same
        # set. Averaging the ratio instead does not: the ratios come in reciprocal pairs of equal
        # weight, and the arithmetic mean of x and 1/x exceeds one unless x = 1.
        @test by_means[126] ≈ 1.0 atol = 1e-12
        @test by_ratios[126] > 1.0
        @test by_ratios[126] ≈ 1.0 atol = 1e-2

        # Away from symmetry the two orders differ, by Jensen's inequality, but only slightly.
        @test by_means[132] != by_ratios[132]
        @test isapprox(by_means[132], by_ratios[132]; rtol = 0.05)
    end

    @testset "weighted mean" begin
        curve = RatioCurve([126, 127], [1.0, 2.0], [1.0, 1.0], "synthetic")
        value, uncertainty = weighted_mean(curve)
        @test value ≈ 1.5
        @test uncertainty ≈ 1 / sqrt(2)
        # Without uncertainties the mean is unweighted and its error is the standard error.
        plain = RatioCurve([126, 127], [1.0, 2.0], [0.0, 0.0], "synthetic")
        @test first(weighted_mean(plain)) ≈ 1.5
        @test_throws ArgumentError weighted_mean(
            RatioCurve(Int[], Float64[], Float64[], "empty")
        )
    end
end
