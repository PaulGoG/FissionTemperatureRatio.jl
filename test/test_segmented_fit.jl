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
