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
        @test segment_count(fit) == 2
        @test only(fit.breakpoints) == 130
        @test first(evaluate(fit, 126)) ≈ 0.5 atol = 0.01
        @test first(evaluate(fit, 170)) ≈ reference_ratio(170) atol = 0.01
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
        # The weight of the median quoted uncertainty, as the consensus imputes it, not the
        # median quoted weight, which differs where the number of quoted values is even.
        @test w[2] ≈ 1 / median([0.1, 0.2])^2
        @test !(w[2] ≈ median([100.0, 25.0]))
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
        f = 52 / 69
        thin = fit_segments(
            masses,
            ratio,
            spread;
            min_segments = 2,
            max_segments = 2,
            measured = fill(f, length(masses)),
        )
        # The same coefficients; χ² and the sample size counted in measurements.
        @test thin.coefficients ≈ full.coefficients
        @test thin.measured_points ≈ f * length(masses)
        @test thin.wrss ≈ f * full.wrss
        @test_throws DimensionMismatch fit_segments(masses, ratio, spread; measured = [1.0])
        @test_throws ArgumentError fit_segments(
            masses, ratio, spread; measured = fill(1.5, length(masses))
        )
        # The fraction scales the information matrix as it scales χ²: the covariance, apart from
        # its chi-squared factor, is that of the full fit over the fraction.
        scale(fit) = max(1.0, fit.wrss / fit.dof)
        @test thin.covariance ./ scale(thin) ≈ full.covariance ./ scale(full) ./ f
    end

    @testset "invalid arguments are rejected" begin
        @test_throws DimensionMismatch fit_segments([1, 2], [1.0], [1.0])
        @test_throws ArgumentError fit_segments(A_H, noisy, σ; max_segments = 0)
        @test_throws ArgumentError fit_segments(A_H, noisy, σ; min_points_per_segment = 1)
        @test_throws ArgumentError fit_segments([3, 2, 1], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0])
        @test_throws ArgumentError fit_weights([0.1, Inf])
        @test_throws ArgumentError fit_weights([0.1, NaN])
        @test_throws "y must be finite" fit_segments(
            A_H, vcat(NaN, noisy[2:end]), σ; max_segments = 2
        )
        @test_throws "y must be finite" fit_segments(
            A_H, vcat(noisy[1:(end - 1)], Inf), σ; max_segments = 2
        )
        @test_throws "pinned_value must be finite" fit_segments(
            A_H, noisy, σ; max_segments = 2, pinned_value = NaN
        )
        @test_throws "bounds must be an ordered open interval" fit_segments(
            A_H, noisy, σ; max_segments = 2, bounds = (1.0, 0.0)
        )
        @test_throws "bounds must be an ordered open interval" fit_segments(
            A_H, noisy, σ; max_segments = 2, bounds = (0.5, 0.5)
        )
    end

    @testset "the number of breakpoint sets a search examines" begin
        count_sets = FissionTemperatureRatio.breakpoint_set_count

        # Every admissible set, by enumeration of the subsets of the interior abscissae.
        function enumerated(x, breakpoints, min_points, windows)
            interior = unique(x)[2:(end - 1)]
            total = 0
            for mask in 0:(2 ^ length(interior) - 1)
                count_ones(mask) == breakpoints || continue
                ψ = [interior[i] for i in eachindex(interior) if (mask >> (i - 1)) & 1 == 1]
                edges = vcat(first(x) - 1, ψ, last(x))
                all(
                    count(xi -> edges[j] < xi ≤ edges[j + 1], x) ≥ min_points for
                    j in 1:(length(edges) - 1)
                ) || continue
                all(window -> any(in(window), ψ), windows) || continue
                total += 1
            end
            return total
        end

        # On n consecutive integers, c breakpoints divide the points into c + 1 runs of at least
        # m: the compositions of n into c + 1 parts of at least m, binomial(n − (c + 1)m + c, c).
        consecutive = collect(126:174)
        @test [count_sets(consecutive, c; min_points_per_segment = 4) for c in 0:5] ==
            [binomial(49 - 4 * (c + 1) + c, c) for c in 0:5] ==
            [1, 42, 741, 7140, 40920, 142506]
        @test count_sets(consecutive, 12; min_points_per_segment = 4) == 0
        @test count_sets(consecutive, 3; min_points_per_segment = 2) isa BigInt

        # Repeated and missing abscissae, with and without windows.
        x = [1, 2, 2, 3, 5, 6, 7, 7, 8, 10, 11, 12, 13, 15, 16]
        for breakpoints in 0:4, min_points in 2:4
            for windows in
                (UnitRange{Int}[], [5:7], [3:4], [2:3, 10:12], [6:6, 7:7], [20:25], [4:4])
                @test count_sets(
                    x,
                    breakpoints;
                    min_points_per_segment = min_points,
                    required_windows = windows,
                ) == enumerated(x, breakpoints, min_points, windows)
            end
        end

        @test_throws ArgumentError count_sets(consecutive, -1)
        @test_throws ArgumentError count_sets(consecutive, 1; min_points_per_segment = 1)
        @test_throws ArgumentError count_sets([3, 2, 1], 1)
    end

    @testset "a search above the limit is refused before it starts" begin
        # Twenty consecutive points, at least four per segment, up to three segments:
        # 1 + 13 + 45 = 59 sets.
        x = collect(120:139)
        y = [a ≤ 130 ? 0.50 - 0.020 * (a - 120) : 0.30 + 0.008 * (a - 130) for a in x]
        y .+= 0.002 .* iseven.(x)
        σ_x = fill(0.004, length(x))
        at_limit = fit_segments(x, y, σ_x; max_segments = 3, max_breakpoint_sets = 59)
        unlimited = fit_segments(x, y, σ_x; max_segments = 3)
        @test at_limit.breakpoints == unlimited.breakpoints
        @test at_limit.coefficients == unlimited.coefficients
        exception = try
            fit_segments(x, y, σ_x; max_segments = 3, max_breakpoint_sets = 58)
        catch caught
            caught
        end
        @test exception isa ArgumentError
        @test occursin("59", exception.msg)
        @test occursin("max_breakpoint_sets", exception.msg)
        # One order alone: 45 sets with two breakpoints.
        @test fit_segments(
            x, y, σ_x; min_segments = 3, max_segments = 3, max_breakpoint_sets = 45
        ) isa SegmentedFit
        @test_throws ArgumentError fit_segments(
            x, y, σ_x; min_segments = 3, max_segments = 3, max_breakpoint_sets = 44
        )
        @test_throws ArgumentError fit_segments(x, y, σ_x; max_breakpoint_sets = 0)
    end

    @testset "pivots span the fitted range" begin
        fit = fit_segments(A_H, noisy, σ; max_segments = 4)
        points = pivots(fit)
        @test first(points)[1] == 126
        @test last(points)[1] == 170
        @test length(points) == segment_count(fit) + 1
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
            @test segment_count(fit) == k
            @test length(fit.selection) == 1
            @test first(only(fit.selection)) == k
        end
        @test segment_count(fit_segments(A_H, value, σ; max_segments = 4)) == 2

        @test_throws ArgumentError fit_segments(A_H, value, σ; min_segments = 0)
        @test_throws ArgumentError fit_segments(
            A_H, value, σ; min_segments = 5, max_segments = 4
        )
    end
end
