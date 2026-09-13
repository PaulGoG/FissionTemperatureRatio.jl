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

    @testset "parsimony biases the choice towards fewer segments" begin
        A_H = collect(120:159)
        value = [a ≤ 130 ? 0.50 - 0.02 * (a - 120) : 0.30 + 0.008 * (a - 130) for a in A_H]
        value .+= 0.004 .* sin.(A_H)
        σ = fill(0.004, length(A_H))

        orders = [
            segments(fit_segments(A_H, value, σ; max_segments = 6, parsimony = λ)) for
            λ in (1.0, 2.0, 4.0, 8.0)
        ]
        # Never increasing with the multiplier. It does not reduce the order here, and should not:
        # the kink in this data is genuine and well resolved, so the likelihood it buys outweighs
        # any plausible penalty. The multiplier bites where the evidence is marginal — on the
        # 252-Cf measurement of Göök it moves the choice from five segments to three.
        @test issorted(orders; rev = true)

        # At one it is the criterion as published, so the default cannot have moved.
        @test segments(fit_segments(A_H, value, σ; max_segments = 6, parsimony = 1.0)) ==
            segments(fit_segments(A_H, value, σ; max_segments = 6))

        # The penalty is what changes, so the criterion reported must change with it even where
        # the chosen order does not.
        plain = fit_segments(A_H, value, σ; max_segments = 6)
        penalised = fit_segments(A_H, value, σ; max_segments = 6, parsimony = 4.0)
        segments(penalised) == segments(plain) && @test penalised.bic > plain.bic

        @test_throws ArgumentError fit_segments(A_H, value, σ; parsimony = 0.0)
        @test_throws ArgumentError fit_segments(A_H, value, σ; parsimony = -1.0)
    end
end
