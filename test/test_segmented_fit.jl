using StableRNGs
using LinearAlgebra: LinearAlgebra

@testset "segmented fit" begin
    rng = StableRNG(20260912)
    A_H = collect(126:170)
    truth = reference_ratio.(A_H)
    noisy = truth .+ 0.004 .* randn(rng, length(A_H))
    σ = fill(0.004, length(A_H))

    @testset "recovers a known breakpoint" begin
        fit = fit_segments(A_H, noisy, σ; max_segments = 4, min_points_per_segment = 4)
        @test segments(fit) == 2
        @test only(fit.breakpoints) == 130
        @test first(evaluate(fit, 126)) ≈ 0.5 atol = 0.01
        @test first(evaluate(fit, 170)) ≈ reference_ratio(170) atol = 0.01
    end

    @testset "continuity holds at every breakpoint" begin
        fit = fit_segments(A_H, noisy, σ; max_segments = 5)
        for ψ in fit.breakpoints
            # Continuity is structural in the truncated-power basis, so the one-sided limits differ
            # only by the step times the slope.
            left = first(evaluate(fit, ψ - 1e-6))
            right = first(evaluate(fit, ψ + 1e-6))
            @test left ≈ right atol = 1e-6
        end
    end

    @testset "pinning is exact and removes the uncertainty there" begin
        fit = fit_segments(A_H, noisy, σ; max_segments = 4, pinned_value = 0.5)
        value, uncertainty = evaluate(fit, 126)
        @test value ≈ 0.5
        @test uncertainty ≈ 0 atol = 1e-12
    end

    @testset "bounds are respected over the whole range" begin
        # A ratio of this form cannot leave (0, 1); a fit that did would make the temperature
        # ratio relation undefined.
        rising = 0.05 .+ 0.02 .* (A_H .- 126)
        fit = fit_segments(
            A_H, rising, fill(0.01, length(A_H)); max_segments = 3, bounds = (0.0, 1.0)
        )
        for x in A_H
            @test 0 < first(evaluate(fit, x)) < 1
        end
        # Data that cannot be described within the bounds leaves no admissible model at all.
        @test_throws InsufficientDataError fit_segments(
            A_H, 2.0 .* rising, fill(0.01, length(A_H)); max_segments = 3, bounds = (0.0, 1.0)
        )
    end

    @testset "the selection criterion is recorded for every order" begin
        fit = fit_segments(A_H, noisy, σ; max_segments = 4)
        @test length(fit.selection) == 4
        @test issorted(first.(fit.selection))
        @test fit.bic == minimum(last.(fit.selection))
    end

    @testset "required windows force a breakpoint" begin
        fit = fit_segments(
            A_H, noisy, σ; max_segments = 4, required_windows = [UnitRange(145, 150)]
        )
        @test any(in(145:150), fit.breakpoints)
    end

    @testset "repeated abscissae are pooled, not treated as separate positions" begin
        doubled_x = repeat(A_H, 2)
        doubled_y = repeat(noisy, 2)
        order = sortperm(doubled_x)
        fit = fit_segments(
            doubled_x[order],
            doubled_y[order],
            fill(0.004, length(doubled_x));
            max_segments = 3,
            min_points_per_segment = 4,
        )
        @test only(fit.breakpoints) == 130
    end

    @testset "weights" begin
        w, imputed = fit_weights([0.1, missing, 0.2])
        @test imputed == 1
        @test w[1] ≈ 100
        @test w[2] ≈ median([100.0, 25.0])
        uniform, all_absent = fit_weights([missing, missing])
        @test all_absent == 2
        @test all(==(1.0), uniform)
        # Zero denotes an exact value, which has no finite weight; an unquoted uncertainty is
        # missing, never zero.
        @test_throws ArgumentError fit_weights([0.1, 0.0])
    end

    @testset "the covariance is scaled by the reduced chi-squared only upwards" begin
        honest = fit_segments(A_H, noisy, σ; max_segments = 3)
        # Uncertainties quoted ten times too large: the residuals are far below what they
        # predict, and the covariance keeps the quoted scale rather than shrinking to the
        # residuals. Uncertainties quoted ten times too small: the covariance inflates.
        generous = fit_segments(A_H, noisy, 10 .* σ; max_segments = 3)
        stingy = fit_segments(A_H, noisy, σ ./ 10; max_segments = 3)
        @test generous.wrss / generous.dof < 1
        @test stingy.wrss / stingy.dof > 1
        @test last(evaluate(generous, 150)) ≈ 10 * last(evaluate(honest, 150)) rtol = 0.2
        @test last(evaluate(stingy, 150)) ≈ last(evaluate(honest, 150)) rtol = 0.2
        # No uncertainty quoted at all: uniform weights carry no scale, so the noise is
        # estimated from the residuals and the result matches the honest fit.
        unquoted = fit_segments(A_H, noisy, fill(missing, length(A_H)); max_segments = 3)
        @test unquoted.weights_imputed == length(A_H)
        @test last(evaluate(unquoted, 150)) ≈ last(evaluate(honest, 150)) rtol = 0.2
    end

    @testset "the covariance of the fitted values" begin
        fit = fit_segments(A_H, noisy, σ; max_segments = 4)
        Σ = covariance(fit, A_H)
        @test size(Σ) == (length(A_H), length(A_H))
        @test Σ ≈ Σ'
        for (index, x) in enumerate(A_H)
            @test sqrt(Σ[index, index]) ≈ last(evaluate(fit, x)) atol = 1e-12
        end
        # Values on one segment are functions of the same two coefficients, so they are
        # correlated, which the diagonal alone cannot say.
        @test abs(Σ[2, 3]) > 0
        @test all(≥(-1e-12), LinearAlgebra.eigvals(LinearAlgebra.Symmetric(Σ)))
    end

    @testset "interpolated points count for their measurements" begin
        masses = collect(126:150)
        ratio = [a ≤ 130 ? 0.5 - 0.03 * (a - 126) : 0.38 + 0.0125 * (a - 130) for a in masses]
        ratio .+= 0.003 .* iseven.(masses)
        spread = fill(0.004, length(masses))
        full = fit_segments(masses, ratio, spread; min_segments = 2, max_segments = 2)
        @test full.points == full.measured_points == length(masses)
        @test full.dof == length(masses) - (length(full.coefficients) + 1)
        f = 52 / 69
        thin = fit_segments(
            masses,
            ratio,
            spread;
            min_segments = 2,
            max_segments = 2,
            measured = fill(f, length(masses)),
        )
        # The same coefficients; χ², sample size and degrees of freedom in measurements.
        @test thin.coefficients ≈ full.coefficients
        @test thin.measured_points ≈ f * length(masses)
        @test thin.wrss ≈ f * full.wrss
        @test thin.dof ≈ f * length(masses) - (length(thin.coefficients) + 1)
        p = length(thin.coefficients) + 1
        n = thin.measured_points
        @test thin.bic ≈ n * log(thin.wrss / n) + p * log(n)
        @test_throws DimensionMismatch fit_segments(masses, ratio, spread; measured = [1.0])
        @test_throws ArgumentError fit_segments(
            masses, ratio, spread; measured = fill(1.5, length(masses))
        )
        # The fraction scales the information matrix as it scales χ²: the covariance, apart from
        # its chi-squared factor, is that of the full fit over the fraction.
        scale(fit) = max(1.0, fit.wrss / fit.dof)
        @test thin.covariance ./ scale(thin) ≈ full.covariance ./ scale(full) ./ f
    end

    @testset "a pool of one interpolated dataset is that dataset" begin
        masses = collect(126:150)
        ratio = reference_ratio.(masses) .+ 0.003 .* iseven.(masses)
        spread = Union{Missing,Float64}[0.003 + 0.0002 * (a - 126) for a in masses]
        spread[7] = missing
        f = 52 / 69
        fractions = fill(f, length(masses))
        curve = RatioCurve(masses, ratio, spread, "interpolated")
        settings = (max_segments = 4, pinned_value = 0.5, bounds = (0.0, 1.0))
        own = fit_segments(masses, ratio, spread; settings..., measured = fractions)

        pooled, measured, σ_measurement = FissionTemperatureRatio._consensus(
            [curve], SYSTEMATIC_TREND_LABEL, [f]
        )
        # The combined curve states the standard error, fraction included; the fit takes the
        # uncertainty of one measurement and the fraction beside it.
        @test measured == fractions
        @test isequal(σ_measurement, spread)
        @test all(skipmissing(pooled.σ .≈ spread ./ sqrt(f)))
        trend = fit_segments(
            pooled.A_H, pooled.ratio, σ_measurement; settings..., measured = measured
        )
        @test segments(trend) > 1
        @test trend.breakpoints == own.breakpoints
        @test trend.coefficients == own.coefficients
        @test trend.wrss == own.wrss
        @test trend.dof == own.dof
        @test trend.bic == own.bic
        @test trend.selection == own.selection
        @test trend.covariance == own.covariance

        # Fitted to the standard error with the fraction beside it, the fraction would enter
        # twice: χ² of a residual at f²/σ², where a measurement of the dataset enters at f/σ².
        order = (min_segments = segments(own), max_segments = segments(own))
        once = fit_segments(masses, ratio, spread; settings..., order..., measured = fractions)
        twice = fit_segments(
            masses, ratio, pooled.σ; settings..., order..., measured = fractions
        )
        @test twice.wrss ≈ f * once.wrss
    end

    @testset "a pool of integer-mass datasets carries no fraction" begin
        masses = collect(126:150)
        ratio = reference_ratio.(masses) .+ 0.003 .* iseven.(masses)
        one = RatioCurve(masses, ratio, fill(0.004, length(masses)), "one")
        other = RatioCurve(
            masses[3:end], ratio[3:end] .+ 0.01, fill(0.006, length(masses) - 2), "other"
        )
        pooled, measured, σ_measurement = FissionTemperatureRatio._consensus(
            [one, other], SYSTEMATIC_TREND_LABEL, [1.0, 1.0]
        )
        # Exactly one, and exactly the standard error: nothing of such a pool moves.
        @test all(==(1.0), measured)
        @test isequal(σ_measurement, pooled.σ)
        counted = fit_segments(
            pooled.A_H, pooled.ratio, σ_measurement; max_segments = 4, measured = measured
        )
        plain = fit_segments(pooled.A_H, pooled.ratio, pooled.σ; max_segments = 4)
        @test counted.coefficients == plain.coefficients
        @test counted.wrss == plain.wrss
        @test counted.covariance == plain.covariance

        # An interpolated dataset beside a measured one: the weight of each combined point is
        # that of its standard error, however it is split between fraction and uncertainty.
        f = 52 / 69
        pooled, measured, σ_measurement = FissionTemperatureRatio._consensus(
            [one, other], SYSTEMATIC_TREND_LABEL, [1.0, f]
        )
        @test measured[1:2] == [1.0, 1.0]
        @test all(f .< measured[3:end] .< 1)
        @test measured ./ σ_measurement .^ 2 ≈ 1 ./ pooled.σ .^ 2
    end

    @testset "a pool of unquoted interpolated datasets carries its fraction once" begin
        masses = collect(126:150)
        n = length(masses)
        base = reference_ratio.(masses)
        unquoted = fill(missing, n)
        one = RatioCurve(masses, base .+ 0.004 .* iseven.(masses), unquoted, "one")
        other = RatioCurve(masses, base .- 0.003 .+ 0.002 .* (masses .% 3), unquoted, "other")
        function trend(factors)
            pooled, measured, σ_measurement = FissionTemperatureRatio._consensus(
                [one, other], SYSTEMATIC_TREND_LABEL, factors
            )
            fit = fit_segments(
                pooled.A_H,
                pooled.ratio,
                σ_measurement;
                min_segments = 2,
                max_segments = 2,
                measured = measured,
            )
            return (; pooled, measured, σ_measurement, fit)
        end
        full = trend([1.0, 1.0])
        half = trend([0.5, 0.5])
        # The standard error is s/√(Σf): it carries the fraction, as σ/√f does for one value,
        # and the uncertainty of one measurement is the same quantity with the fraction out.
        @test all(==(0.5), half.measured)
        @test half.pooled.ratio == full.pooled.ratio
        @test half.pooled.σ ≈ full.pooled.σ ./ sqrt(0.5)
        @test half.σ_measurement ≈ full.pooled.σ
        # So χ² scales with the fraction, once, and the degrees of freedom count measurements.
        @test half.fit.coefficients ≈ full.fit.coefficients
        @test half.fit.wrss ≈ 0.5 * full.fit.wrss
        @test half.fit.dof ≈ 0.5 * n - (length(half.fit.coefficients) + 1)
        @test half.fit.dof < full.fit.dof

        # Unequal fractions: the weight of every combined point is that of its standard error.
        mixed = trend([1.0, 0.5])
        @test all(≈((1.0 + 0.25) / 1.5), mixed.measured)
        @test mixed.pooled.σ ≈ full.pooled.σ .* sqrt(2 / 1.5)
        @test mixed.measured ./ mixed.σ_measurement .^ 2 ≈ 1 ./ mixed.pooled.σ .^ 2
    end

    @testset "identical unquoted values leave the combined point unquoted" begin
        masses = collect(126:150)
        base = reference_ratio.(masses)
        unquoted = fill(missing, length(masses))
        shifted = base .+ 0.004
        # At one mass the two datasets tabulate the same value: there is no dispersion to
        # estimate an uncertainty from, and a zero would claim an exact value and stop the fit.
        shifted[5] = base[5]
        one = RatioCurve(masses, base, unquoted, "one")
        other = RatioCurve(masses, shifted, unquoted, "other")
        pooled, measured, σ_measurement = FissionTemperatureRatio._consensus(
            [one, other], SYSTEMATIC_TREND_LABEL, [1.0, 1.0]
        )
        @test ismissing(pooled.σ[5]) && ismissing(σ_measurement[5])
        @test count(ismissing, pooled.σ) == 1
        @test pooled.ratio[5] == base[5]
        fit = fit_segments(
            pooled.A_H, pooled.ratio, σ_measurement; max_segments = 3, measured = measured
        )
        @test fit.weights_imputed == 1
    end

    @testset "the between-dataset variance is estimated against the expectation of Q" begin
        # Two values ±d about their mean, σ = 1, each counting for f of a measurement: the
        # weights are f, Q = 2 f d² and its expectation without a variance between the datasets
        # is Σf − Σf²/Σf = f, not k − 1 = 1. With d² = 0.75 and f = 1/2, Q = 0.75 exceeds 0.5,
        # τ² = (0.75 − 0.5)/(1 − 0.5) = 0.5, and the standard error is (Σ f/(1 + τ²))^(-1/2).
        d = sqrt(0.75)
        low = RatioCurve([130], [0.3 - d], [1.0], "low")
        high = RatioCurve([130], [0.3 + d], [1.0], "high")
        interpolated = consensus([low, high]; weights = [0.5, 0.5])
        @test only(interpolated.ratio) ≈ 0.3
        @test only(interpolated.σ) ≈ sqrt(1.5)
        # Measured at integer masses, Q = 1.5 against k − 1 = 1: DerSimonian and Laird's
        # estimate, τ² = 0.5, unchanged.
        measured = consensus([low, high])
        @test only(measured.σ) ≈ sqrt(1.5 / 2)
    end

    @testset "invalid arguments are rejected" begin
        @test_throws DimensionMismatch fit_segments([1, 2], [1.0], [1.0])
        @test_throws ArgumentError fit_segments(A_H, noisy, σ; max_segments = 0)
        @test_throws ArgumentError fit_segments(A_H, noisy, σ; min_points_per_segment = 1)
        @test_throws ArgumentError fit_segments([3, 2, 1], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0])
    end

    @testset "pivots span the fitted range" begin
        fit = fit_segments(A_H, noisy, σ; max_segments = 4)
        points = pivots(fit)
        @test first(points)[1] == 126
        @test last(points)[1] == 170
        @test length(points) == segments(fit) + 1
    end

    @testset "a single order can be fitted rather than selected" begin
        A_H = collect(120:139)
        value = [a ≤ 130 ? 0.50 - 0.02 * (a - 120) : 0.30 + 0.008 * (a - 130) for a in A_H]
        value .+= 0.002 .* iseven.(A_H)
        σ = fill(0.004, length(A_H))

        # Setting the two bounds equal fits exactly that many segments, which is how one order is
        # inspected on its own; the criterion still prefers two here.
        for k in 1:4
            fit = fit_segments(A_H, value, σ; min_segments = k, max_segments = k)
            @test segments(fit) == k
            @test length(fit.selection) == 1
            @test first(only(fit.selection)) == k
        end
        @test segments(fit_segments(A_H, value, σ; max_segments = 4)) == 2

        @test_throws ArgumentError fit_segments(A_H, value, σ; min_segments = 0)
        @test_throws ArgumentError fit_segments(
            A_H, value, σ; min_segments = 5, max_segments = 4
        )
    end

    @testset "a minimum span rejects segments shorter than it" begin
        A_H = collect(120:159)
        value = [a ≤ 130 ? 0.50 - 0.02 * (a - 120) : 0.30 + 0.008 * (a - 130) for a in A_H]
        value .+= 0.004 .* sin.(A_H)
        σ = fill(0.004, length(A_H))

        # With two points per segment the point count alone admits a segment one mass unit
        # long; the span guard is what keeps every segment at least six.
        fit = fit_segments(
            A_H, value, σ; max_segments = 6, min_points_per_segment = 2, min_segment_span = 6
        )
        pivot_abscissae = [first(A_H); fit.breakpoints; last(A_H)]
        @test all(≥(6), diff(pivot_abscissae))
        # At the shipped point count on consecutive mass numbers a span of three coincides with
        # the point guard, so it changes nothing.
        @test fit_segments(A_H, value, σ; max_segments = 4, min_segment_span = 3).breakpoints ==
            fit_segments(A_H, value, σ; max_segments = 4).breakpoints
        # A span the data cannot honour leaves no admissible model.
        @test_throws InsufficientDataError fit_segments(
            A_H, value, σ; min_segments = 2, max_segments = 2, min_segment_span = 30
        )
        @test_throws ArgumentError fit_segments(A_H, value, σ; min_segment_span = 0)
    end
end
