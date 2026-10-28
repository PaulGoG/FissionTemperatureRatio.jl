# The combination of datasets into the consensus curve: the pooling fractions of interpolated
# datasets, the dispersion of unquoted pools and the estimate of the between-dataset variance.

using StableRNGs: StableRNG

@testset "consensus" begin
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

    @testset "the dispersion of an unquoted pool is taken with the weights of its mean" begin
        # Values 0, 1, 3 with f = 1, 1/2, 1/2: the mean is (0 + 0.5 + 1.5)/2 = 1,
        # Σ f (r − r̄)² = 1 + 0 + 2 = 3, its expectation in units of s² is Σf − Σf²/Σf = 1.25,
        # so s² = 2.4 and the standard error is √(s²/Σf) = √1.2. The unweighted sample variance
        # about the unweighted mean, 7/3, would give another figure.
        curves = [RatioCurve([130], [r], [missing], "set $(r)") for r in (0.0, 1.0, 3.0)]
        combined = consensus(curves; weights = [1.0, 0.5, 0.5])
        @test only(combined.ratio) ≈ 1.0
        @test only(combined.σ) ≈ sqrt(1.2)
        # Equal fractions: the sample variance, and s/√(Σf).
        equal = consensus(curves; weights = [0.5, 0.5, 0.5])
        @test only(equal.σ) ≈ sqrt((7 / 3) / 1.5)
    end

    @testset "the estimate of the between-dataset variance is unbiased under pooling factors" begin
        # Six values of variance σ² + τ² about one mean, combined with the factors f: the
        # estimate of τ² averages to τ² over many draws. Against k − 1 in place of the
        # expectation of Q it would average about σ²(1 − f)/f lower, 0.84 here. With equal σ
        # the estimate is recovered from the standard error, Σ f/(σ² + τ²) = 1/σ_r̄².
        rng = StableRNG(20261002)
        quoted, τ² = 1.0, 4.0
        factors = [1.0, 0.5, 0.5, 0.4, 0.6, 0.5]
        draws = 40_000
        estimates = Vector{Float64}(undef, draws)
        uncertainties = Union{Missing,Float64}[quoted for _ in factors]
        for draw in 1:draws
            values = sqrt(quoted^2 + τ²) .* randn(rng, length(factors))
            _, spread, _, _ = FissionTemperatureRatio._combine(values, uncertainties, factors)
            estimates[draw] = spread^2 * sum(factors) - quoted^2
        end
        # The truncation at zero adds about 0.02 at these values; the Monte Carlo standard error
        # of the mean is about 0.02.
        @test sum(estimates) / draws ≈ τ² atol = 0.1
    end

    @testset "identical unquoted values of unequal fractions leave the point unquoted" begin
        # The weighted mean of identical values differs from them by rounding where the
        # fractions are unequal; that is no dispersion, and no uncertainty is quoted from it.
        unquoted = Union{Missing,Float64}[missing, missing]
        for factors in ([0.3, 0.7], [0.5, 0.5], [1.0, 0.417])
            value, spread, _, single = FissionTemperatureRatio._combine(
                [0.1, 0.1], unquoted, factors
            )
            @test value ≈ 0.1
            @test ismissing(spread)
            @test ismissing(single)
        end
        # Values that differ in their seventh significant digit do have a dispersion.
        _, spread, _, _ = FissionTemperatureRatio._combine(
            [0.1, 0.1000001], unquoted, [0.3, 0.7]
        )
        @test spread > 1e-9
        # The shares of the combination are the normalized fractions.
        @test last(FissionTemperatureRatio._combine([0.1, 0.2], unquoted, [0.3, 0.7])) ≈
            [0.3, 0.7]
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
end
