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
- `outside_physical_range`: points whose ratio falls outside `(0, 1)`, where the transformation to
  the temperature ratio is undefined.
- `symmetry_departure`: `|r_ν(A₀/2) - 1/2|`, exactly zero for a faithful measurement, or `missing`
  where the set has no symmetric split.
- `coverage`: the fraction of the mass numbers of the fragmentation range at which the dataset
  provides a complete pair. A dataset below the configured floor is diagnosed but offers no
  segmented curve: with pairs at few mass numbers the breakpoint search cannot place the minimum
  where the data do not reach, and the curve it returns asserts structure between measurements.
- `complement_sum`, `complement_spread`: mean and standard deviation of `ν(A) + ν(A₀-A)` over the
  pairs. A per-fragment dataset sums to about the total multiplicity with small spread; a large
  spread indicates a normalization or a quantity problem.
- `without_uncertainties`: points quoting no uncertainty, which are given the median weight.
"""
struct DatasetDiagnostics
    label::String
    points::Int
    pairs::Int
    first_pair::Union{Int,Missing}
    last_pair::Union{Int,Missing}
    coverage::Float64
    outside_physical_range::Int
    symmetry_departure::Union{Float64,Missing}
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
        ismissing(complement) && continue
        push!(sums, data.ν[index] + complement[1])
    end

    symmetric = if iseven(A₀)
        index = findfirst(==(A₀ ÷ 2), ratio.A_H)
        index === nothing ? missing : abs(ratio.ratio[index] - 0.5)
    else
        missing
    end

    return DatasetDiagnostics(
        data.label,
        length(data),
        length(ratio),
        isempty(ratio) ? missing : first(ratio.A_H),
        isempty(ratio) ? missing : last(ratio.A_H),
        count(in(A_H_range), ratio.A_H) / length(A_H_range),
        count(v -> v ≤ 0 || v ≥ 1, ratio.ratio),
        symmetric,
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
w = 1 / (σ² + τ²),    r̄ = Σ w r / Σ w,    σ_r̄ = (Σ w)^(-1/2).
```

Concatenating the datasets and weighting by the quoted uncertainties alone would let the dataset
with the smallest quoted uncertainties decide the result, and would weight a dataset by its
number of points. Because the datasets disagree by ten to twenty times their quoted
uncertainties, `τ²` dominates, the weights become nearly equal, and the uncertainty of the
combination reflects the disagreement.

Where a mass number has one measurement only, that value and its uncertainty pass through: there
is no dispersion to estimate. Where none of the values at a mass number carries a quoted
uncertainty, the mean is taken and the standard error of the values supplies the uncertainty,
`s/√k` with `s` their standard deviation about it; identical values leave nothing to estimate it
from, and the combined point then quotes no uncertainty. Points quoting no uncertainty alongside
points that do take the median of the quoted uncertainties at that mass number;
[`fit_weights`](@ref) imputes the median weight instead, which differs from the median
uncertainty where the quoted count is even.

`weights`, one factor per curve, scales the weight of every value of that curve, in the
fixed-effect stage and in the combination alike: `w = f/(σ² + τ²)`. A dataset whose masses were
interpolated onto the integers has rows that share their bracketing points and so are not
independent; [`pooling_weight`](@ref) gives it the factor raw points over rows. A curve alone at a
mass number passes through with its uncertainty divided by `√f`, and values quoting no
uncertainty are averaged with the weights `f` to the standard error `s/√(Σf)`, `s²` being
`Σ f (r − r̄)²/(Σf − Σf²/Σf)`.

The standard error written here reads `f` as the share of a measurement a value amounts to. The
estimate of `τ²` reads it as a weight on a value of variance `σ² + τ²`; for the mean under that
reading the written standard error is conservative, by `1/√f` where the fractions are equal.

With such factors the fixed-effect weights are `f/σ²` while a value keeps the variance `σ²`, and
the statistic `Q = Σ (f/σ²)(r − r̄)²` has the expectation `Σf − Σ(f²/σ²)/Σ(f/σ²)` where the
datasets do not differ, not `k − 1`. `τ²` is estimated against that expectation, the method of
moments for general weights of DerSimonian and Kacker, Contemp. Clin. Trials **28**, 105 (2007),
doi:10.1016/j.cct.2006.04.004; it reduces to DerSimonian and Laird's where every factor is one.

The uncertainty of the curve returned is the standard error of each combined value, factor
included. A fit to it must not weight a point by its measured fraction a second time; the
systematic trend is fitted to the same values at the uncertainty of one measurement, with the
fraction as [`fit_segments`](@ref)'s `measured`, which applies it once.
"""
function consensus(
    curves::Vector{RatioCurve};
    label::AbstractString = SYSTEMATIC_TREND_LABEL,
    weights::AbstractVector{<:Real} = ones(length(curves)),
)
    return first(_consensus(curves, label, weights))
end

# The combined curve, the measured fraction of each point, the uncertainty of each point as one
# measurement, and per point the members combined there, as indices into `curves`, with their
# shares of the combination.
function _consensus(curves::Vector{RatioCurve}, label::AbstractString, weights::AbstractVector)
    length(weights) == length(curves) || throw(
        DimensionMismatch("one weight per curve: $(length(weights)) for $(length(curves))")
    )
    for (curve, weight) in zip(curves, weights)
        weight isa Real ||
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
            push!(factors, weight isa Real ? weight : weight[index])
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

    return RatioCurve(A_H, ratio, σ, String(label)), measured, σ_measurement, members, shares
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
    # One value: its own uncertainty is that of one measurement, exactly.
    k == 1 && return (
        values[1], uncertainties[1] / sqrt(factors[1]), factors[1], uncertainties[1], [1.0]
    )

    quoted = .!ismissing.(uncertainties)
    if !any(quoted)
        # No value quotes an uncertainty: the weighted mean, with the dispersion s of the values
        # about it standing for the uncertainty of one of them.
        total = sum(factors)
        μ = sum(factors .* values) / total
        fraction = sum(factors .^ 2) / total
        dispersion = sqrt(sum(factors .* (values .- μ) .^ 2) / (total - fraction))
        # Identical values leave no dispersion to estimate from.
        share = factors ./ total
        dispersion > UNRESOLVED_DISPERSION * maximum(abs, values) ||
            return (μ, missing, fraction, missing, share)
        spread = dispersion / sqrt(total)
        return (μ, spread, fraction, spread * sqrt(fraction), share)
    end

    σ = Float64[coalesce(u, NaN) for u in uncertainties]
    σ[.!quoted] .= median(view(σ, quoted))

    # Fixed-effect combination first, since the between-dataset variance is estimated from its
    # residuals.
    w₀ = factors ./ σ .^ 2
    total = sum(w₀)
    fixed = sum(w₀ .* values) / total
    Q = sum(w₀ .* (values .- fixed) .^ 2)
    # E[Q] with no between-dataset variance, for the weights f/σ²: Σf − Σ w₀f/Σ w₀ (DerSimonian
    # and Kacker); k − 1 only where every f is one.
    expected = sum(factors) - sum(w₀ .* factors) / total
    # The estimator is truncated at zero: a Q below its expectation means the datasets agree
    # better than their uncertainties suggest, not that the variance between them is negative.
    τ² = max(0.0, (Q - expected) / (total - sum(w₀ .^ 2) / total))

    w = factors ./ (σ .^ 2 .+ τ²)
    spread = 1 / sqrt(sum(w))
    fraction = sum(w .* factors) / sum(w)
    return (sum(w .* values) / sum(w), spread, fraction, spread * sqrt(fraction), w ./ sum(w))
end
