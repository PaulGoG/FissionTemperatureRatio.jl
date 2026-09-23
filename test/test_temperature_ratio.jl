
@testset "temperature ratio" begin
    @testset "equal sharing gives equal temperatures" begin
        # r_ν = 1/2 with R_a = 1 means the pair shares its excitation energy equally between two
        # fragments of equal level density parameter, so their temperatures coincide.
        curve = RatioCurve([126], [0.5], [0.0], "synthetic")
        R_T = temperature_ratio(curve, Dict(126 => 1.0))
        @test only(R_T.ratio) ≈ 1.0
    end

    @testset "relation is the inverse of the excitation energy partition" begin
        R_a = 1.3
        for r in (0.2, 0.35, 0.5, 0.7)
            curve = RatioCurve([130], [r], [0.0], "synthetic")
            R_T = only(temperature_ratio(curve, Dict(130 => R_a)).ratio)
            # E*_H/TXE = 1/(1 + R_a R_T²) must return the ratio it was built from.
            @test 1 / (1 + R_a * R_T^2) ≈ r
        end
    end

    @testset "uncertainty propagation" begin
        r, σ_r, R_a = 0.3, 0.02, 1.2
        curve = RatioCurve([130], [r], [σ_r], "synthetic")
        result = temperature_ratio(curve, Dict(130 => R_a))
        R_T = only(result.ratio)
        @test only(result.σ) ≈ σ_r / (2 * R_T * R_a * r^2)
    end

    @testset "singular and absent points are omitted" begin
        curve = RatioCurve([126, 130], [0.5, 0.4], [0.0, 0.0], "synthetic")
        @test length(temperature_ratio(curve, Dict(126 => 1.0))) == 1
    end

    DATA_AVAILABLE && @testset "averaging over the charge distribution" begin
        A₀, Z₀ = 252, 98
        domain = fragmentation_domain(A₀, Z₀, 126:132, 5, FLAT_CHARGES)
        model = BSFG_MODEL
        by_means = level_density_ratio(RatioOfMeans(), model, A₀, Z₀, domain)
        by_ratios = level_density_ratio(MeanOfRatios(), model, A₀, Z₀, domain)

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

    @testset "an unquoted uncertainty stays unquoted" begin
        curve = RatioCurve([130, 132], [0.4, 0.45], [missing, 0.01], "synthetic")
        result = temperature_ratio(curve, Dict(130 => 1.2, 132 => 1.1))
        @test ismissing(result.σ[1])
        @test result.σ[2] > 0
    end

    @testset "a segmented curve carries the covariance of its temperature ratio" begin
        A_H = collect(126:150)
        value = [a ≤ 130 ? 0.5 - 0.03 * (a - 126) : 0.38 + 0.0125 * (a - 130) for a in A_H]
        value .+= 0.002 .* iseven.(A_H)
        fit = fit_segments(
            A_H, value, fill(0.004, length(A_H)); max_segments = 3, pinned_value = 0.5
        )
        R_a = Dict(a => (252 - a) / a for a in A_H)
        curve = SegmentedCurve("synthetic", fit, R_a)

        @test curve.R_T.A_H == A_H
        @test size(curve.R_T_covariance) == (length(A_H), length(A_H))
        # The tabulated pointwise uncertainty is the diagonal of the propagated covariance.
        for index in eachindex(A_H)
            @test sqrt(curve.R_T_covariance[index, index]) ≈ curve.R_T.σ[index] atol = 1e-12
        end
        # Exact at the pin: zero uncertainty there, and R_T = 1 by the synthetic R_a.
        @test curve.R_T.σ[1] == 0
        @test curve.R_T.ratio[1] ≈ 1

        yields = MassYield(
            A_H, exp.(-((A_H .- 140) ./ 5) .^ 2), fill(missing, length(A_H)), "y", ""
        )
        value_cov, σ_cov = total_average(curve, yields)
        value_ind, σ_ind = total_average(curve.R_T, yields)
        # Same value; the covariance form carries the correlation the independent form drops.
        @test value_cov ≈ value_ind
        @test σ_cov > σ_ind
        record = TotalAverage(curve, yields)
        @test record.value == value_cov
        @test record.uncertainty == σ_cov
        @test record.uncertainty_independent_points == σ_ind

        mean, σ = range_mean(curve)
        @test mean ≈ sum(curve.R_T.ratio) / length(A_H)
        @test σ > 0
        # A mass number the level density ratio does not cover is dropped, covariance included.
        partial = SegmentedCurve("partial", fit, Dict(a => R_a[a] for a in 126:140))
        @test partial.R_T.A_H == collect(126:140)
        @test size(partial.R_T_covariance) == (15, 15)
        @test_throws ArgumentError SegmentedCurve("none", fit, Dict(200 => 1.0))
    end
end
