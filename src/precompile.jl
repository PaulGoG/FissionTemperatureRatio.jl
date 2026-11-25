# The calls a first extraction compiles, run once on synthetic curves when the package is
# precompiled: the combination of two datasets, the correlogram of their deviations and the
# correlation matrix formed from it, and the segmented fit under that matrix with its evaluation.
# Nothing is read or written.

@setup_workload begin
    A_H = collect(126:150)
    ratio = [a ≤ 131 ? 0.50 - 0.02 * (a - 126) : 0.40 + 0.01 * (a - 131) for a in A_H]
    σ = Union{Missing,Float64}[0.004 for _ in A_H]
    @compile_workload begin
        whole = RatioCurve(A_H, ratio .+ 0.003 .* iseven.(A_H), σ, "whole")
        interpolated = RatioCurve(A_H, ratio .- 0.002 .* iseven.(A_H), σ, "interpolated")
        curves = [whole, interpolated]
        fractions = [ones(length(A_H)), fill(0.7, length(A_H))]
        pooled = _consensus(curves, SYSTEMATIC_TREND_LABEL, fractions)
        correlogram = deviation_correlogram(curves, pooled.curve; weights = mean.(fractions))
        autocorrelation_decay(correlogram; lags = DEFAULT_AUTOCORRELATION_LAGS)
        R = pooled_correlation(pooled.curve.A_H, pooled.members, pooled.shares, 0.5)
        fit = fit_segments(
            pooled.curve.A_H,
            pooled.curve.ratio,
            pooled.σ_measurement;
            max_segments = 3,
            measured = pooled.measured,
            pinned_value = 0.5,
            required_windows = [129:133],
            bounds = (0.0, 1.0),
            correlation = R,
            max_breakpoint_sets = DEFAULT_MAX_BREAKPOINT_SETS,
        )
        evaluate(fit, 135)
        covariance(fit, A_H)
        pivots(fit)
        unresolved_minimum(fit)
    end
end
