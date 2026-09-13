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
    outside_physical_range::Int
    symmetry_departure::Union{Float64,Missing}
    complement_sum::Union{Float64,Missing}
    complement_spread::Union{Float64,Missing}
    without_uncertainties::Int
end

"""
    diagnose(data, ratio, A₀) -> DatasetDiagnostics

Describe one multiplicity dataset and the ratio extracted from it.

Computes only what can be judged without reference to the other datasets; see
[`DatasetDiagnostics`](@ref) for why that restriction is deliberate.
"""
function diagnose(data::Multiplicity, ratio::RatioCurve, A₀::Integer)
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
        count(v -> v ≤ 0 || v ≥ 1, ratio.ratio),
        symmetric,
        isempty(sums) ? missing : mean(sums),
        length(sums) < 2 ? missing : std(sums),
        count(≤(0), data.σν),
    )
end

"""
    consensus(curves; label) -> RatioCurve

Combine several measurements of the same ratio into one curve, mass number by mass number.

At each mass number the available values are combined by inverse-variance weighting with an
additional between-dataset variance `τ²`, estimated from their dispersion after DerSimonian and
Laird,
Control. Clin. Trials **7**, 177 (1986):

```
w = 1 / (σ² + τ²),    r̄ = Σ w r / Σ w,    σ_r̄ = (Σ w)^(-1/2).
```

This is what the pooled curve requires rather than a refinement of it. Concatenating the datasets
and weighting by the quoted uncertainties alone hands the result to whichever author quoted the
smallest ones, and counts a dataset with many points more heavily than one with few, neither of
which is a statement about the measurements. Because the datasets here disagree by ten to twenty
times their quoted uncertainties, `τ²` dominates, the weights become nearly equal, and the
uncertainty of the combination reflects the disagreement instead of hiding it.

Where a mass number has one measurement only, that value and its uncertainty pass through: there
is no dispersion to estimate. Where none of the values at a mass number carries an uncertainty,
the unweighted mean is taken and the standard error of the values supplies the uncertainty.
Points quoting no uncertainty alongside points that do are given the median of the latter, as
[`fit_weights`](@ref) does.
"""
function consensus(curves::Vector{RatioCurve}; label::AbstractString = TREND_LABEL)
    masses = sort!(unique!(reduce(vcat, (curve.A_H for curve in curves); init = Int[])))
    A_H = Int[]
    ratio = Float64[]
    σ = Float64[]

    for mass in masses
        values = Float64[]
        uncertainties = Float64[]
        for curve in curves
            index = findfirst(==(mass), curve.A_H)
            index === nothing && continue
            push!(values, curve.ratio[index])
            push!(uncertainties, curve.σ[index])
        end
        isempty(values) && continue

        combined, spread = _combine(values, uncertainties)
        push!(A_H, mass)
        push!(ratio, combined)
        push!(σ, spread)
    end

    return RatioCurve(A_H, ratio, σ, String(label))
end

function _combine(values::Vector{Float64}, uncertainties::Vector{Float64})
    k = length(values)
    k == 1 && return (values[1], uncertainties[1])

    positive = uncertainties .> 0
    if !any(positive)
        μ = mean(values)
        return (μ, std(values) / sqrt(k))
    end

    σ = copy(uncertainties)
    σ[.!positive] .= median(view(uncertainties, positive))

    # Fixed-effect combination first, since the between-dataset variance is estimated from its
    # residuals.
    w₀ = 1 ./ σ .^ 2
    total = sum(w₀)
    fixed = sum(w₀ .* values) / total
    Q = sum(w₀ .* (values .- fixed) .^ 2)
    # The estimator is truncated at zero: a Q below its expectation means the datasets agree
    # better than their uncertainties suggest, not that the variance between them is negative.
    τ² = max(0.0, (Q - (k - 1)) / (total - sum(w₀ .^ 2) / total))

    w = 1 ./ (σ .^ 2 .+ τ²)
    return (sum(w .* values) / sum(w), 1 / sqrt(sum(w)))
end
