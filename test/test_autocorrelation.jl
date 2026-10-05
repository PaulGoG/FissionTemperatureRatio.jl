using StableRNGs
using LinearAlgebra: LinearAlgebra
using Statistics: mean, var

# The reference numbers of this file come from a dense evaluation of the formulae, independent of
# the package: Σ = D R D with D = W^(-1/2), M = I − X(XᵀWX)⁻¹XᵀW, the covariance
# (XᵀWX)⁻¹XᵀWΣWX(XᵀWX)⁻¹ and the expectation tr(MᵀWMΣ) of the residual sum of squares.

@testset "autocorrelated errors" begin
    unquoted(n) = Vector{Union{Missing,Float64}}(fill(missing, n))
    # The correlation matrix of an AR(1) error along the abscissa.
    ar1(x, ρ) = [ρ^abs(a - b) for a in x, b in x]

    @testset "the correlogram of the deviations" begin
        combined = RatioCurve(collect(1:5), fill(0.5, 5), unquoted(5), "combined")
        # Deviations less their mean, in units of 0.01: −0.8, 0.2, 1.2, 0.2, −0.8, and the
        # opposite for the second dataset.
        rising = [1.0, 2.0, 3.0, 2.0, 1.0]
        above = RatioCurve(collect(1:5), 0.5 .+ 0.01 .* rising, unquoted(5), "above")
        below = RatioCurve(collect(1:5), 0.5 .- 0.01 .* rising, unquoted(5), "below")
        ρ = deviation_correlogram([above, below], combined)
        @test length(ρ) == CORRELOGRAM_LAGS == 8
        # Per lag: the products over the pairs present, over the root of the squares of their
        # leading and of their trailing members.
        @test ρ[1] ≈ 0.16 / 2.16
        @test ρ[2] ≈ -1.88 / 2.12
        @test ρ[3] ≈ -0.32 / 0.68
        @test ρ[4] ≈ 1.0
        @test all(ismissing, ρ[5:8])
        @test length(deviation_correlogram([above, below], combined; lags = 3)) == 3
        @test_throws ArgumentError deviation_correlogram([above, below], combined; lags = 0)
        @test_throws DimensionMismatch deviation_correlogram(
            [above, below], combined; weights = [1.0]
        )

        # A dataset on a two-unit grid has no pair one or three mass units apart: it says
        # nothing about the odd lags and leaves them where they were.
        alternate = RatioCurve(
            [1, 3, 5], 0.5 .+ 0.01 .* [1.0, -2.0, 1.0], unquoted(3), "two-unit grid"
        )
        with_grid = deviation_correlogram([above, below, alternate], combined)
        @test with_grid[1] ≈ ρ[1]
        @test with_grid[3] ≈ ρ[3]
        @test with_grid[2] ≈ (2 * -1.88 - 4.0) / (2 * 2.12 + 5.0)
        # A dataset enters with the factor it is pooled with.
        weighted = deviation_correlogram(
            [above, below, alternate], combined; weights = [1.0, 1.0, 0.5]
        )
        @test weighted[2] ≈ (2 * -1.88 - 2.0) / (2 * 2.12 + 2.5)
        @test weighted[1] ≈ ρ[1]

        # A mass number one dataset holds alone is left out: its deviation there is zero by
        # construction, whatever its error.
        longer = RatioCurve(collect(1:7), vcat(above.ratio, [0.9, 0.1]), unquoted(7), "longer")
        wider = RatioCurve(collect(1:7), vcat(fill(0.5, 5), [0.9, 0.1]), unquoted(7), "wider")
        @test deviation_correlogram([longer, below], wider)[1:4] ≈ ρ[1:4]
        @test all(ismissing, deviation_correlogram([above], combined))

        # A constant offset of a dataset is removed with its mean.
        shifted = RatioCurve(above.A_H, above.ratio .+ 0.1, unquoted(5), "shifted")
        @test deviation_correlogram([shifted, below], combined)[1:4] ≈ ρ[1:4]

        # Nothing to correlate: datasets identical to the combined curve, or sharing one mass.
        same = RatioCurve(collect(1:5), fill(0.5, 5), unquoted(5), "same")
        single = RatioCurve([3], [0.6], unquoted(1), "single")
        @test all(ismissing, deviation_correlogram([same, same], combined))
        @test all(ismissing, deviation_correlogram([single, above], combined))
        @test all(ismissing, deviation_correlogram(RatioCurve[], combined))
    end

    @testset "the correlation of a combined curve's errors" begin
        # Two datasets at every mass, a third at the odd ones, equal shares where they meet.
        masses = collect(1:6)
        members = [isodd(A) ? [1, 2, 3] : [1, 2] for A in masses]
        shares = [fill(1 / length(m), length(m)) for m in members]
        ρ = 0.6
        R = pooled_correlation(masses, members, shares, ρ)
        @test R ≈ R'
        @test all(R[m, m] == 1 for m in eachindex(masses))
        # Neighbours share two datasets of the three and two held there: 2/√6 of ρ.
        @test R[1, 2] ≈ ρ * 2 / sqrt(6)
        # Next-neighbours hold the same datasets with the same shares: ρ² undiminished.
        @test R[1, 3] ≈ ρ^2
        @test R[2, 4] ≈ ρ^2
        @test R[1, 4] ≈ ρ^3 * 2 / sqrt(6)
        @test minimum(LinearAlgebra.eigvals(LinearAlgebra.Symmetric(R))) > 0
        # The same datasets everywhere: the AR(1) matrix itself.
        everywhere = pooled_correlation(masses, fill([1, 2], 6), fill([0.3, 0.7], 6), ρ)
        @test everywhere ≈ ar1(masses, ρ)
        # No dataset in common, or no correlation along the axis: independent points.
        apart = pooled_correlation([1, 2], [[1], [2]], [[1.0], [1.0]], ρ)
        @test apart == [1.0 0.0; 0.0 1.0]
        @test pooled_correlation(masses, members, shares, 0.0) == LinearAlgebra.I(6)
        # A gap in the mass numbers counts in the power.
        @test pooled_correlation([1, 4], [[1], [1]], [[1.0], [1.0]], ρ)[1, 2] ≈ ρ^3
        @test_throws ArgumentError pooled_correlation(masses, members, shares, 1.0)
        @test_throws ArgumentError pooled_correlation([1], [[1, 2]], [[0.5, 0.4]], ρ)
        @test_throws DimensionMismatch pooled_correlation([1, 2], [[1]], [[1.0]], ρ)
    end

    @testset "the decay fitted to the first lags" begin
        geometric = [0.5^k for k in 1:8]
        @test autocorrelation_decay(geometric) ≈ 0.5 atol = 1e-6
        @test autocorrelation_decay(geometric; lags = 1) ≈ 0.5 atol = 1e-6
        @test autocorrelation_decay(geometric; lags = 8) ≈ 0.5 atol = 1e-6
        # A second lag above the first: the first lag alone gives 0.4427, the first four 0.6894.
        staggered = [0.4427, 0.5290, 0.2999, 0.3880]
        @test autocorrelation_decay(staggered; lags = 1) ≈ 0.4427 atol = 1e-6
        @test autocorrelation_decay(staggered; lags = 2) ≈ 0.620890 atol = 2e-6
        @test autocorrelation_decay(staggered; lags = 4) ≈ 0.689449 atol = 2e-6
        # More lags asked for than the correlogram holds: those it holds.
        @test autocorrelation_decay(staggered; lags = 6) ≈ 0.689449 atol = 2e-6
        # A lag no dataset spans is left out of the sum.
        @test autocorrelation_decay([missing, 0.25]; lags = 2) ≈ 0.5 atol = 1e-6
        # Bounded below by zero and above by 0.999.
        @test autocorrelation_decay([-0.9, 0.8, -0.7]) == 0.0
        @test autocorrelation_decay([-0.3]) == 0.0
        @test autocorrelation_decay(fill(1.0, 4)) ≈ 0.999 atol = 1e-6
        @test autocorrelation_decay([missing, missing, missing, missing]) == 0.0
        @test autocorrelation_decay(Union{Missing,Float64}[]) == 0.0
        @test_throws ArgumentError autocorrelation_decay(geometric; lags = 0)
    end

    A_H = collect(120:139)
    r = [a ≤ 130 ? 0.50 - 0.020 * (a - 120) : 0.30 + 0.008 * (a - 130) for a in A_H]
    r .+= 0.002 .* iseven.(A_H)
    σ = fill(0.004, length(A_H))

    @testset "independent points are the default, to the bit" begin
        default = fit_segments(A_H, r, σ; max_segments = 3)
        explicit = fit_segments(A_H, r, σ; max_segments = 3, correlation = nothing)
        @test !default.correlated
        @test default.expected_wrss == default.dof
        for field in fieldnames(SegmentedFit)
            @test isequal(getfield(default, field), getfield(explicit, field))
        end
        # The identity as a correlation matrix: the same covariance, by the other route.
        identity = fit_segments(A_H, r, σ; max_segments = 3, correlation = ar1(A_H, 0.0))
        @test identity.correlated
        @test identity.covariance ≈ default.covariance rtol = 1e-12
        @test identity.expected_wrss ≈ default.dof rtol = 1e-12
    end

    @testset "the kernel changes the covariance and nothing else" begin
        independent = fit_segments(A_H, r, σ; max_segments = 3)
        correlated = fit_segments(A_H, r, σ; max_segments = 3, correlation = ar1(A_H, 0.6))
        for field in (:breakpoints, :coefficients, :wrss, :dof, :bic, :selection, :points)
            @test getfield(correlated, field) == getfield(independent, field)
        end
        @test correlated.correlated
        @test correlated.breakpoints == [130]
        # Dense reference: tr(MᵀWMΣ) = 10.979589328 of n − p = 17, on dof = 16.
        @test correlated.wrss ≈ 1.24038896427 rtol = 1e-9
        @test correlated.expected_wrss ≈ 10.3337311323 rtol = 1e-9
        @test correlated.expected_wrss < correlated.dof
        @test last(evaluate(correlated, 125)) ≈ 0.00210195917509 rtol = 1e-8
        @test last(evaluate(correlated, 135)) ≈ 0.00222876174861 rtol = 1e-8
        @test correlated.covariance[1, 1] ≈ 1.27940352109e-05 rtol = 1e-8
        @test correlated.covariance[end, end] ≈ 9.54053566821e-07 rtol = 1e-8
        @test correlated.covariance ≈ correlated.covariance'
        @test last(evaluate(independent, 125)) ≈ 0.00113940163782 rtol = 1e-8
    end

    @testset "the scale refers to the expectation of wrss under the kernel" begin
        # Uncertainties a quarter of the above: wrss is sixteen times larger, 19.85, which is 1.24
        # times dof = 16 and 1.92 times its expectation under the kernel.
        tight = fill(0.001, length(A_H))
        independent = fit_segments(A_H, r, tight; max_segments = 3)
        correlated = fit_segments(A_H, r, tight; max_segments = 3, correlation = ar1(A_H, 0.6))
        wrss = 16 * 1.24038896427
        @test correlated.wrss ≈ wrss rtol = 1e-9
        @test correlated.expected_wrss ≈ 10.3337311323 rtol = 1e-9
        @test last(evaluate(correlated, 125)) ≈
            0.00210195917509 / 4 * sqrt(wrss / 10.3337311323) rtol = 1e-8
        @test last(evaluate(independent, 125)) ≈ 0.00113940163782 / 4 * sqrt(wrss / 16) rtol =
            1e-8
    end

    @testset "with the pin and with measured fractions" begin
        pinned = fit_segments(
            A_H, r, σ; max_segments = 3, pinned_value = 0.5, correlation = ar1(A_H, 0.6)
        )
        @test pinned.breakpoints == [130]
        @test pinned.dof == 17
        @test pinned.expected_wrss ≈ 12.9192454707 rtol = 1e-9
        @test last(evaluate(pinned, 125)) ≈ 0.00154109502308 rtol = 1e-8
        @test last(evaluate(pinned, 135)) ≈ 0.00221520291928 rtol = 1e-8
        @test last(evaluate(pinned, 120)) ≈ 0 atol = 1e-12

        # A fraction common to every point: the degrees of freedom count measurements, the share
        # of them a correlated error leaves in the residuals is the same.
        halved = fit_segments(
            A_H, r, σ; max_segments = 3, measured = fill(0.5, 20), correlation = ar1(A_H, 0.6)
        )
        @test halved.dof == 6
        @test halved.expected_wrss ≈ 3.87514917461 rtol = 1e-9
        @test last(evaluate(halved, 125)) ≈ 0.00297261917297 rtol = 1e-8
    end

    @testset "a long series: the variance of the mean grows by (1 + ρ)/(1 − ρ)" begin
        x = collect(1:400)
        line = 0.2 .+ 0.001 .* x
        quoted = fill(0.01, 400)
        independent = fit_segments(x, line, quoted; max_segments = 1)
        correlated = fit_segments(x, line, quoted; max_segments = 1, correlation = ar1(x, 0.6))
        # The fit of a line passes through the mean at the centroid; 1 + 2 Σ (1 − k/n) ρᵏ.
        factor = (last(evaluate(correlated, 200.5)) / last(evaluate(independent, 200.5)))^2
        @test factor ≈ 3.98125 rtol = 1e-8
        @test factor ≈ (1 + 0.6) / (1 - 0.6) rtol = 0.01
    end

    @testset "simulated AR(1) errors: residuals and scatter as predicted" begin
        rng = StableRNG(20261005)
        n, ρ, s = 60, 0.6, 0.01
        x = collect(1:n)
        quoted = fill(s, n)
        centroid = (n + 1) / 2
        R = ar1(x, ρ)
        residual_sums = Float64[]
        centroid_values = Float64[]
        expected = NaN
        for _ in 1:3000
            ε = Vector{Float64}(undef, n)
            ε[1] = s * randn(rng)
            for t in 2:n
                ε[t] = ρ * ε[t - 1] + s * sqrt(1 - ρ^2) * randn(rng)
            end
            fit = fit_segments(
                x, 0.2 .+ 0.005 .* x .+ ε, quoted; max_segments = 1, correlation = R
            )
            push!(residual_sums, fit.wrss)
            push!(centroid_values, first(evaluate(fit, centroid)))
            expected = fit.expected_wrss
            if length(residual_sums) == 1
                scale = max(1.0, fit.wrss / fit.expected_wrss)
                predicted = s^2 / n * (1 + 2 * sum((1 - k / n) * ρ^k for k in 1:(n - 1)))
                @test last(evaluate(fit, centroid))^2 / scale ≈ predicted rtol = 1e-8
            end
        end
        # One segment, unpinned, every point one measurement: dof = n − 2 and the expectation is
        # tr(PR) itself.
        @test expected < n - 2
        @test mean(residual_sums) ≈ expected rtol = 0.03
        @test var(centroid_values) ≈
            s^2 / n * (1 + 2 * sum((1 - k / n) * ρ^k for k in 1:(n - 1))) rtol = 0.1
        # Against the degrees of freedom the residuals would look too small.
        @test mean(residual_sums) < n - 2 - 3
    end

    @testset "simulated datasets on different grids: the combined curve as predicted" begin
        # Three datasets with AR(1) errors of one variance, the third on the odd masses only,
        # combined by their mean at each mass: the errors of the combined points are those
        # `pooled_correlation` describes, and a line fitted to them scatters accordingly.
        rng = StableRNG(20261006)
        n, ρ, s = 40, 0.7, 0.01
        x = collect(1:n)
        members = [isodd(A) ? [1, 2, 3] : [1, 2] for A in x]
        held = length.(members)
        shares = [fill(1 / k, k) for k in held]
        R = pooled_correlation(x, members, shares, ρ)
        quoted = s ./ sqrt.(held)
        function ar1_errors()
            ε = Vector{Float64}(undef, n)
            ε[1] = s * randn(rng)
            for t in 2:n
                ε[t] = ρ * ε[t - 1] + s * sqrt(1 - ρ^2) * randn(rng)
            end
            return ε
        end
        centroid_values = Float64[]
        residual_sums = Float64[]
        weights = held ./ s^2
        centroid = sum(weights .* x) / sum(weights)
        predicted = NaN
        expected = NaN
        for _ in 1:3000
            ε = [ar1_errors() for _ in 1:3]
            combined = [sum(ε[i][A] for i in members[A]) / held[A] for A in x]
            fit = fit_segments(
                x, 0.2 .+ 0.005 .* x .+ combined, quoted; max_segments = 1, correlation = R
            )
            push!(centroid_values, first(evaluate(fit, centroid)))
            push!(residual_sums, fit.wrss)
            scale = max(1.0, fit.wrss / fit.expected_wrss)
            predicted = last(evaluate(fit, centroid))^2 / scale
            expected = fit.expected_wrss
        end
        @test var(centroid_values) ≈ predicted rtol = 0.1
        @test mean(residual_sums) ≈ expected rtol = 0.03
        # With every pair of points held at ρ^|ΔA| the scatter would be overstated, and with
        # independent points understated.
        stationary = fit_segments(
            x, 0.2 .+ 0.005 .* x, quoted; max_segments = 1, correlation = ar1(x, ρ)
        )
        independent = fit_segments(x, 0.2 .+ 0.005 .* x, quoted; max_segments = 1)
        @test last(evaluate(stationary, centroid))^2 > 1.05 * predicted
        @test last(evaluate(independent, centroid))^2 < 0.5 * predicted
    end

    @testset "arguments" begin
        @test_throws DimensionMismatch fit_segments(A_H, r, σ; correlation = ar1(1:5, 0.5))
        skewed = ar1(A_H, 0.5)
        skewed[1, 2] = 0.1
        @test_throws ArgumentError fit_segments(A_H, r, σ; correlation = skewed)
        scaled = 2 .* ar1(A_H, 0.5)
        @test_throws ArgumentError fit_segments(A_H, r, σ; correlation = scaled)
        beyond = ar1(A_H, 0.5)
        beyond[1, 2] = beyond[2, 1] = 1.5
        @test_throws ArgumentError fit_segments(A_H, r, σ; correlation = beyond)
    end

    @testset "a required window is not met by a curve without breakpoints" begin
        @test_throws InsufficientDataError fit_segments(
            A_H, r, σ; max_segments = 1, required_windows = [128:132]
        )
        # Out of reach of every candidate: no model satisfies it, the single segment included.
        @test_throws InsufficientDataError fit_segments(
            A_H, r, σ; max_segments = 3, required_windows = [150:155]
        )
        @test fit_segments(A_H, r, σ; max_segments = 3, required_windows = [128:132]).breakpoints ==
            [130]
        @test segments(fit_segments(A_H, r, σ; max_segments = 1)) == 1
    end

    @testset "a curve that resolves no minimum says so" begin
        turning = fit_segments(A_H, r, σ; max_segments = 3)
        @test unresolved_minimum(turning) === nothing

        ripple = 0.0005 .* iseven.(A_H)
        rising = fit_segments(A_H, 0.2 .+ 0.01 .* (A_H .- 120) .+ ripple, σ; max_segments = 1)
        @test occursin("rises from A_H = 120", unresolved_minimum(rising))
        falling = fit_segments(A_H, 0.6 .- 0.01 .* (A_H .- 120) .+ ripple, σ; max_segments = 1)
        @test occursin("falls to A_H = 139", unresolved_minimum(falling))

        # Two segments that both rise: a kink, not a minimum.
        kinked = [a ≤ 130 ? 0.20 + 0.002 * (a - 120) : 0.22 + 0.02 * (a - 130) for a in A_H]
        fit = fit_segments(A_H, kinked .+ ripple, σ; max_segments = 3)
        @test segments(fit) ≥ 2
        @test occursin("rises from A_H = 120", unresolved_minimum(fit))
    end
end
