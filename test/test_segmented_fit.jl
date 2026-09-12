using StableRNGs

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
        @test_throws ArgumentError fit_segments(
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
        w, imputed = fit_weights([0.1, 0.0, 0.2])
        @test imputed == 1
        @test w[1] ≈ 100
        @test w[2] ≈ median([100.0, 25.0])
        uniform, all_absent = fit_weights([0.0, 0.0])
        @test all_absent == 2
        @test all(==(1.0), uniform)
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
end
