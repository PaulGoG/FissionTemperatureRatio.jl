# One segmented description of the multiplicity ratio, carried through to the temperature ratio
# with the covariance the transformation gives it.

"""
    ExtractedCurve

One extracted temperature ratio: the segmented fit of the multiplicity ratio, the temperature
ratio derived from it, and the diagnostics of that curve.

A run produces several: one per experimental dataset, and one following the systematic behaviour
of the ratio. They are alternatives, not an ensemble to be averaged — a prompt emission code takes
one of them as input, and which one describes reality is settled downstream, by comparing the
multiplicity distributions and yields that code produces against experiment. The curve a
consuming code reads is the tabulated `R_T`, through the manifest a run writes; it is not this
object.

Construct with `ExtractedCurve(label, fit, fitted_to, averaging, model, domain)`, which evaluates
the fit at every mass number of its range, inverts the relation between `R_T` and `E*_H/TXE` there
as `averaging` prescribes, and propagates the coefficient covariance through the exact slope of
that inverse.

# Fields

- `label`: the dataset the curve came from, or `SYSTEMATIC_TREND_LABEL`.
- `kind`: `"dataset"` or `"systematic_trend"`, as the manifest declares it.
- `fit`: the [`SegmentedFit`](@ref) of `r_ν`.
- `r_ν`, `R_T`: the segmented multiplicity ratio and the temperature ratio from it, tabulated at
  every mass number of the fitted range where the relation is defined.
- `R_T_covariance`: covariance of `R_T.ratio`. Its diagonal is the square of `R_T.σ`; the
  off-diagonal elements are what a yield-weighted average of the curve has to carry, since a
  curve with a handful of coefficients cannot have independent errors at every mass number.
- `pairs`: the points the fit was made to: the dataset's complete fragment pairs, or the mass
  numbers of the combined curve.
- `coverage`: the fraction of the mass numbers of the fragmentation range those points reach.
"""
struct ExtractedCurve
    label::String
    kind::String
    fit::SegmentedFit
    r_ν::RatioCurve
    R_T::RatioCurve
    R_T_covariance::Matrix{Float64}
    pairs::Int
    coverage::Float64

    function ExtractedCurve(
        label::AbstractString,
        fit::SegmentedFit,
        fitted_to::RatioCurve,
        averaging::RatioAveraging,
        model::LevelDensityModel,
        domain::FragmentationDomain,
    )
        r_ν = evaluate(fit, fit.x₀:fit.x_max; label = label)
        R_T = temperature_ratio(averaging, model, domain, r_ν)
        isempty(R_T) && throw(
            ArgumentError(
                "the temperature ratio relation is defined at no mass number of the range \
                 $(fit.x₀):$(fit.x_max) of the curve $(repr(label))"
            ),
        )
        # Covariance of r_ν at the mass numbers R_T kept, then the inverse's slope at each.
        Σ = covariance(fit, R_T.A_H)
        # Defined wherever the inverse is, and R_T keeps only those mass numbers.
        slope = Float64[
            something(temperature_ratio_slope(averaging, model, domain, A, R)) for
            (A, R) in zip(R_T.A_H, R_T.ratio)
        ]
        range = domain.heavy_masses
        kind = label == SYSTEMATIC_TREND_LABEL ? "systematic_trend" : "dataset"
        coverage = count(in(range), unique(fitted_to.A_H)) / length(range)
        return new(
            String(label),
            kind,
            fit,
            r_ν,
            R_T,
            Matrix(Symmetric(slope .* Σ .* slope')),
            length(fitted_to),
            coverage,
        )
    end
end

# The variance of a linear functional `wᵀ R_T` of the curve.
function _functional_variance(curve::ExtractedCurve, w::AbstractVector{<:Real})
    return max(dot(w, curve.R_T_covariance, w), 0.0)
end

"""
    range_mean(curve) -> Tuple{Float64,Float64}

The unweighted mean of the temperature ratio of a segmented curve over its mass range, with the
uncertainty propagated through the curve's covariance.

This is *not* the total average the literature quotes: it weights every mass number equally and
so is dominated by the far-asymmetric tail, where the yield is negligible. It is the only summary
available when no yield distribution is given.
"""
function range_mean(curve::ExtractedCurve)
    n = length(curve.R_T)
    w = fill(1 / n, n)
    return (dot(w, curve.R_T.ratio), sqrt(_functional_variance(curve, w)))
end

"""
    total_average(curve::ExtractedCurve, yields) -> Tuple{Float64,Float64}

The total average `⟨R_T⟩` of a segmented curve over a fragment mass yield distribution, with the
fit covariance propagated.

The value is that of [`total_average`](@ref) applied to the curve's tabulated temperature ratio.
The uncertainty differs: the tabulated points are functions of a few fitted coefficients, so the
ratio term is `wᵀ C w` with `w` the normalized yield weights and `C` the curve's covariance,
rather than the sum of squares the independent-points form takes. The yield term is unchanged.
For the curves of the four shipped systems the covariance-propagated uncertainty is two to five
times the independent-points one, which is the published approximation and remains available as
`total_average(curve.R_T, yields)`.

Throws as [`total_average`](@ref) does when the curve and the distribution share no mass number.
"""
function FissionFragmentsDomain.total_average(curve::ExtractedCurve, yields::MassYield)
    R_T = curve.R_T
    weight = zeros(length(R_T))
    σ_Y = zeros(length(R_T))
    for (index, mass) in enumerate(R_T.A_H)
        entry = mass_yield(yields, mass)
        entry === nothing && continue
        weight[index] = entry[1]
        σ_Y[index] = coalesce(entry[2], 0.0)
    end

    total = sum(weight)
    any(mass_yield(yields, mass) !== nothing for mass in R_T.A_H) ||
        throw(ArgumentError("the ratio curve \"$(curve.label)\" and the yield distribution \
             \"$(yields.label)\" share no mass number"))
    total > 0 || throw(
        ArgumentError(
            "the yield distribution \"$(yields.label)\" sums to $(total) over the mass \
             numbers it shares with \"$(curve.label)\""
        ),
    )

    w = weight ./ total
    mean = dot(w, R_T.ratio)
    variance = _functional_variance(curve, w)
    for index in eachindex(w)
        variance += ((R_T.ratio[index] - mean) / total)^2 * σ_Y[index]^2
    end
    return (mean, sqrt(variance))
end

"""
    TotalAverage

The total average of one segmented curve over one yield distribution: the value, its uncertainty
with the fit covariance propagated, the uncertainty of the independent-points form, which is the
approximation the published tables use and is kept for comparison with them, and the fraction of
the distribution's yield over the fragmentation range that the curve takes in,
[`yield_fraction`](@ref).

    TotalAverage(curve, yields, heavy_masses)

The total average of `curve` over `yields`, `heavy_masses` being the heavy mass numbers of the
fragmentation range.
"""
struct TotalAverage
    value::Float64
    uncertainty::Float64
    uncertainty_independent_points::Float64
    yield_fraction::Float64
end

function TotalAverage(
    curve::ExtractedCurve, yields::MassYield, heavy_masses::AbstractUnitRange{<:Integer}
)
    value, uncertainty = total_average(curve, yields)
    return TotalAverage(
        value,
        uncertainty,
        last(total_average(curve.R_T, yields)),
        yield_fraction(curve.R_T, yields, heavy_masses),
    )
end
