# Combining several measurements of the same ratio, and the diagnostics that describe each of them.

"""
    DatasetDiagnostics

What can be said about one experimental dataset without judging its values against the others.

The distinction matters here. Datasets of one fissioning system disagree far beyond their quoted
uncertainties — for 252-Cf the spread between datasets at a given mass number runs to ten or
twenty times the median quoted uncertainty — so a criterion of the form "reject what is more than
a few sigma from the consensus" rejects everything, and a reduced chi-squared ranks how generously
an author quoted errors rather than how good the measurement is. These diagnostics are therefore
structural: coverage, and departures from identities the ratio satisfies by construction.

# Fields

- `label`: the dataset.
- `points`: multiplicity values read.
- `pairs`: complete fragment pairs within the heavy-mass range, which is what a fit actually has
  to work with.
- `first_pair`, `last_pair`: the heavy-mass span those pairs cover.
- `outside_physical_range`: complete fragment pairs within the heavy-mass range whose ratio falls
  outside `(0, 1)`, where the transformation to the temperature ratio is undefined, and which the
  ratio therefore leaves out: a fragment of the pair at zero multiplicity.
- `coverage`: the fraction of the mass numbers of the fragmentation range at which the dataset
  provides a complete pair. A dataset below the configured floor is diagnosed but offers no
  segmented curve: with pairs at few mass numbers the breakpoint search cannot place the minimum
  where the data do not reach, and the curve it returns asserts structure between measurements.
- `complement_sum`, `complement_spread`: mean and standard deviation of `ν(A) + ν(A₀-A)` over the
  pairs. A per-fragment dataset sums to about the total multiplicity with small spread; a large
  spread indicates a normalization or a quantity problem.
- `without_uncertainties`: points quoting no uncertainty, which are given the median of the
  quoted uncertainties.
"""
struct DatasetDiagnostics
    label::String
    points::Int
    pairs::Int
    first_pair::Union{Int,Missing}
    last_pair::Union{Int,Missing}
    coverage::Float64
    outside_physical_range::Int
    complement_sum::Union{Float64,Missing}
    complement_spread::Union{Float64,Missing}
    without_uncertainties::Int
end

"""
    diagnose(data, ratio, A₀, A_H_range) -> DatasetDiagnostics

Describe one multiplicity dataset and the ratio extracted from it over the fragmentation range
`A_H_range`.

Computes only what can be judged without reference to the other datasets; see
[`DatasetDiagnostics`](@ref) for why that restriction is deliberate.
"""
function diagnose(data::Multiplicity, ratio::RatioCurve, A₀::Integer, A_H_range::UnitRange{Int})
    sums = Float64[]
    for (index, A) in enumerate(data.A)
        2 * A ≤ A₀ && continue
        complement = multiplicity(data, A₀ - A)
        complement === nothing && continue
        push!(sums, data.ν[index] + complement[1])
    end

    # The pairs the ratio leaves out: both fragments measured, the ratio outside (0, 1).
    outside = 0
    for A_H in A_H_range
        heavy = multiplicity(data, A_H)
        light = multiplicity(data, A₀ - A_H)
        (heavy === nothing || light === nothing) && continue
        total = heavy[1] + light[1]
        (total > 0 && 0 < heavy[1] / total < 1) || (outside += 1)
    end

    return DatasetDiagnostics(
        data.label,
        length(data),
        length(ratio),
        isempty(ratio) ? missing : first(ratio.A_H),
        isempty(ratio) ? missing : last(ratio.A_H),
        count(in(A_H_range), ratio.A_H) / length(A_H_range),
        outside,
        isempty(sums) ? missing : mean(sums),
        length(sums) < 2 ? missing : std(sums),
        count(ismissing, data.σν),
    )
end

"""
    consensus(curves; label, weights) -> RatioCurve

Combine several measurements of the same ratio into one curve, mass number by mass number.

At each mass number the available values are combined by inverse-variance weighting with an
additional between-dataset variance `τ²`, estimated from their dispersion after DerSimonian and
Laird, Control. Clin. Trials **7**, 177 (1986), doi:10.1016/0197-2456(86)90046-2:

```
v = σ²/f,    w = 1/(v + τ²),    r̄ = Σ w r / Σ w,    σ_r̄ = (Σ w)^(-1/2).
```

Concatenating the datasets and weighting by the quoted uncertainties alone would let the dataset
with the smallest quoted uncertainties decide the result, and would weight a dataset by its
number of points. Because the datasets disagree by ten to twenty times their quoted
uncertainties, `τ²` dominates, the weights become nearly equal, and the uncertainty of the
combination reflects the disagreement.

`weights` gives each curve its measured fraction `f`, the share of a measurement each of its
values amounts to, one by default. A dataset whose masses were interpolated onto the integers
has rows that share their bracketing points and so are not independent; its
[`RetrievalRecord`](@ref) gives it its measured values over its rows. Such a value carries the
variance of one measurement inflated by the share it amounts to, `v = σ²/f`, while `τ²`, by which
a dataset's curve departs from the others as a whole, is not diluted by the interpolation. The
fixed-effect weights `1/v` are then inverse variances, the statistic
`Q = Σ (r − r̄)²/v` has the expectation `k − 1` where the datasets do not differ, and `τ²` is
DerSimonian and Laird's moment estimate, truncated at zero.

Where a mass number has one measurement only, that value passes through with the uncertainty
`σ/√f`: there is no dispersion to estimate. Where none of the values at a mass number carries a
quoted uncertainty, each is taken to have the variance `s²/f` of its share of a measurement of
variance `s²`, the mean is weighted by `f`, and `s² = Σ f (r − r̄)²/(k − 1)` supplies the
uncertainty `s/√(Σf)`, the standard error of the mean where the fractions are equal; identical
values leave nothing to estimate it from, and the combined point then quotes no uncertainty.
Values quoting no uncertainty beside values that do take the median of the quoted uncertainties
at that mass number, as [`fit_weights`](@ref) does.

The uncertainty of the curve returned is the standard error of each combined value, measured
fraction included. A fit to it must not weight a point by its measured fraction a second time;
the systematic trend is fitted to the same values at the uncertainty of one measurement, with the
fraction as [`fit_segments`](@ref)'s `measured`, which applies it once.
"""
function consensus(
    curves::Vector{RatioCurve};
    label::AbstractString = SYSTEMATIC_TREND_LABEL,
    weights::AbstractVector{<:Real} = ones(length(curves)),
)
    length(weights) == length(curves) || throw(
        DimensionMismatch("one weight per curve: $(length(weights)) for $(length(curves))")
    )
    per_point = [fill(Float64(w), length(curve)) for (curve, w) in zip(curves, weights)]
    return _consensus(curves, label, per_point).curve
end

# A combination of several curves, mass number by mass number: the combined curve, the measured
# fraction of each point, the uncertainty of each point as one measurement, and per point the
# members combined there, as indices into the curves, with their shares of the combination.
struct Combination
    curve::RatioCurve
    measured::Vector{Float64}
    σ_measurement::Vector{Union{Missing,Float64}}
    members::Vector{Vector{Int}}
    shares::Vector{Vector{Float64}}
end

# The combination of `curves`, `weights` holding the measured fraction of every point of each.
function _consensus(
    curves::Vector{RatioCurve},
    label::AbstractString,
    weights::AbstractVector{<:AbstractVector{<:Real}},
)
    length(weights) == length(curves) || throw(
        DimensionMismatch("one weight per curve: $(length(weights)) for $(length(curves))")
    )
    for (curve, weight) in zip(curves, weights)
        length(weight) == length(curve) ||
            throw(DimensionMismatch("one pooling weight per point of $(repr(curve.label)): \
                     $(length(weight)) for $(length(curve))"))
        all(f -> 0 < f <= 1, weight) ||
            throw(ArgumentError("a pooling weight lies in (0, 1], got $(weight)"))
    end
    masses = sort!(unique!(reduce(vcat, (curve.A_H for curve in curves); init = Int[])))
    A_H = Int[]
    ratio = Float64[]
    σ = Union{Missing,Float64}[]
    measured = Float64[]
    σ_measurement = Union{Missing,Float64}[]
    members = Vector{Int}[]
    shares = Vector{Float64}[]

    for mass in masses
        values = Float64[]
        uncertainties = Union{Missing,Float64}[]
        factors = Float64[]
        present = Int[]
        for (position, (curve, weight)) in enumerate(zip(curves, weights))
            index = findfirst(==(mass), curve.A_H)
            index === nothing && continue
            push!(values, curve.ratio[index])
            push!(uncertainties, curve.σ[index])
            push!(factors, weight[index])
            push!(present, position)
        end
        isempty(values) && continue

        combined, spread, fraction, single, share = _combine(values, uncertainties, factors)
        push!(A_H, mass)
        push!(ratio, combined)
        push!(σ, spread)
        push!(measured, fraction)
        push!(σ_measurement, single)
        push!(members, present)
        push!(shares, share)
    end

    return Combination(
        RatioCurve(A_H, ratio, σ, String(label)), measured, σ_measurement, members, shares
    )
end

# The dispersion of unquoted values, relative to the largest of them, below which they are taken
# as identical: far above rounding, far below the last digit a tabulation carries.
const UNRESOLVED_DISPERSION = 1e-10

function _combine(
    values::Vector{Float64},
    uncertainties::Vector{Union{Missing,Float64}},
    factors::Vector{Float64} = ones(length(values)),
)
    k = length(values)
    # One value: its own uncertainty is that of one measurement, exactly, and its variance that
    # of the share of a measurement it amounts to.
    k == 1 && return (
        values[1], uncertainties[1] / sqrt(factors[1]), factors[1], uncertainties[1], [1.0]
    )

    quoted = .!ismissing.(uncertainties)
    if !any(quoted)
        # No value quotes an uncertainty: each has the variance s²/f of its share of one
        # measurement, the mean is weighted by f, and E[Σ f (r − r̄)²] = s² (k − 1) gives s².
        total = sum(factors)
        μ = sum(factors .* values) / total
        fraction = sum(factors .^ 2) / total
        share = factors ./ total
        dispersion = sqrt(sum(factors .* (values .- μ) .^ 2) / (k - 1))
        # Identical values leave no dispersion to estimate from.
        dispersion > UNRESOLVED_DISPERSION * maximum(abs, values) ||
            return (μ, missing, fraction, missing, share)
        spread = dispersion / sqrt(total)
        return (μ, spread, fraction, spread * sqrt(fraction), share)
    end

    σ = Float64[coalesce(u, NaN) for u in uncertainties]
    σ[.!quoted] .= median(view(σ, quoted))
    # The variance of each value as the share f of a measurement of variance σ²; the fixed-effect
    # weights are its inverse, so Q has the expectation k − 1 without a between-dataset variance.
    v = σ .^ 2 ./ factors
    a = 1 ./ v
    total = sum(a)
    fixed = sum(a .* values) / total
    Q = sum(a .* (values .- fixed) .^ 2)
    # The estimator is truncated at zero: a Q below its expectation means the datasets agree
    # better than their uncertainties suggest, not that the variance between them is negative.
    τ² = max(0.0, (Q - (k - 1)) / (total - sum(a .^ 2) / total))

    w = 1 ./ (v .+ τ²)
    spread = 1 / sqrt(sum(w))
    fraction = sum(w .* factors) / sum(w)
    return (sum(w .* values) / sum(w), spread, fraction, spread * sqrt(fraction), w ./ sum(w))
end
