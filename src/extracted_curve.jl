# One segmented description of the multiplicity ratio, carried through to the temperature ratio
# with the covariance the transformation gives it.

# The delete-one jackknife covariance of the values of `refits` at `masses`, the Gram matrix of
# their scaled deviations: Σ[a, b] = Σᵢ dᵢ(a) dᵢ(b), with dᵢ(a) = √((m − 1)/m) (θᵢ(a) − θ̄(a))
# over the m refits that hold the mass a. A refit that does not hold a mass contributes nothing
# there, and neither does a mass fewer than two refits hold.
function _jackknife_covariance(
    masses::AbstractVector{<:Integer}, refits::AbstractVector{RatioCurve}
)
    deviation = zeros(length(refits), length(masses))
    for (column, mass) in enumerate(masses)
        held = Int[]
        θ = Float64[]
        for (row, refit) in enumerate(refits)
            index = findfirst(==(mass), refit.A_H)
            index === nothing && continue
            push!(held, row)
            push!(θ, refit.ratio[index])
        end
        m = length(θ)
        m < 2 && continue
        θ̄ = mean(θ)
        scale = sqrt((m - 1) / m)
        for (row, value) in zip(held, θ)
            deviation[row, column] = scale * (value - θ̄)
        end
    end
    return Matrix(Symmetric(deviation' * deviation))
end

"""
    ExtractedCurve

One extracted temperature ratio: the segmented fit of the multiplicity ratio, the temperature
ratio derived from it, and the diagnostics of that curve.

A run produces several: one per experimental dataset, and one following the systematic behaviour
of the ratio. They are alternatives; a consuming code takes one, and which one describes reality
is settled downstream, by comparing the multiplicity distributions and yields that code produces
against experiment. The curve a
consuming code reads is the tabulated `R_T`, through the manifest a run writes; it is not this
object.

Construct with `ExtractedCurve(label, fit, fitted_to, averaging, model, domain)`, which evaluates
the fit at every mass number of its range, inverts the relation between `R_T` and `E*_H/TXE` there
as `averaging` prescribes, and propagates the coefficient covariance through the exact slope of
that inverse.

`ExtractedCurve(curve, refits)` gives `curve` with the delete-one jackknife over `refits`, which
are the same curve fitted with each pooled dataset left out in turn: `r_ν.σ` and `R_T.σ` become
the larger of the propagated uncertainty and the jackknife at each mass number, and
`R_T_jackknife_covariance` holds the jackknife covariance of `R_T`. The values, the fit and
`R_T_covariance` are those of `curve`. It throws an `ArgumentError` for fewer than two refits.
The method page gives the estimator and what it was calibrated against.

# Fields

- `label`: the dataset the curve came from, or `SYSTEMATIC_TREND_LABEL`, which makes it the
  systematic trend; the manifest declares the kind of every curve so.
- `fit`: the [`SegmentedFit`](@ref) of `r_ν`; `fit.points` is the number of points it was made
  to, the dataset's complete fragment pairs or the mass numbers of the combined curve.
- `r_ν`, `R_T`: the segmented multiplicity ratio and the temperature ratio from it, tabulated at
  every mass number of the fitted range where the relation is defined.
- `R_T_covariance`: the covariance of `R_T.ratio` propagated from the fit. Its off-diagonal
  elements are what a yield-weighted average of the curve has to carry, since a curve with a
  handful of coefficients cannot have independent errors at every mass number. Its diagonal is
  the square of `R_T.σ`, except for a curve that carries a jackknife, where `R_T.σ` is the
  larger of the two.
- `coverage`: the fraction of the mass numbers of the fragmentation range those points reach.
- `R_T_jackknife_covariance`: the delete-one jackknife covariance of `R_T.ratio` over the
  leave-one-out refits, for the systematic trend of a run that made at least two; `nothing` for
  every other curve. Where it is present, the uncertainty of `r_ν` and of `R_T` at each mass
  number, and that of every average of the curve, is the larger of the propagated one and the
  jackknife of the same quantity.
"""
struct ExtractedCurve
    label::String
    fit::SegmentedFit
    r_ν::RatioCurve
    R_T::RatioCurve
    R_T_covariance::Matrix{Float64}
    coverage::Float64
    R_T_jackknife_covariance::Union{Nothing,Matrix{Float64}}

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
        coverage = count(in(range), unique(fitted_to.A_H)) / length(range)
        return new(
            String(label),
            fit,
            r_ν,
            R_T,
            Matrix(Symmetric(slope .* Σ .* slope')),
            coverage,
            nothing,
        )
    end

    function ExtractedCurve(curve::ExtractedCurve, refits::AbstractVector{ExtractedCurve})
        length(refits) < 2 && throw(
            ArgumentError(
                "the jackknife of a curve takes at least two refits, got $(length(refits))"
            ),
        )
        Σ_r = _jackknife_covariance(curve.r_ν.A_H, RatioCurve[refit.r_ν for refit in refits])
        Σ_J = _jackknife_covariance(curve.R_T.A_H, RatioCurve[refit.R_T for refit in refits])
        # A pinned value is exact, and a refit that reaches the pin holds it as well: what the
        # refits differ by there is rounding, and the curve keeps a zero uncertainty at its pin.
        if curve.fit.pinned_value !== nothing
            for (Σ, masses) in ((Σ_r, curve.r_ν.A_H), (Σ_J, curve.R_T.A_H))
                pin = findfirst(==(curve.fit.x₀), masses)
                pin === nothing && continue
                Σ[pin, :] .= 0.0
                Σ[:, pin] .= 0.0
            end
        end
        σ_r = sqrt.(diag(Σ_r))
        σ_R = sqrt.(diag(Σ_J))
        r_ν = RatioCurve(
            copy(curve.r_ν.A_H),
            copy(curve.r_ν.ratio),
            Vector{Union{Missing,Float64}}(max.(curve.r_ν.σ, σ_r)),
            curve.r_ν.label,
        )
        R_T = RatioCurve(
            copy(curve.R_T.A_H),
            copy(curve.R_T.ratio),
            Vector{Union{Missing,Float64}}(max.(curve.R_T.σ, σ_R)),
            curve.R_T.label,
        )
        return new(curve.label, curve.fit, r_ν, R_T, curve.R_T_covariance, curve.coverage, Σ_J)
    end
end

# The kind of a curve as the manifest declares it: `"systematic_trend"` for the trend, `"dataset"`
# for any other.
_is_trend(curve::ExtractedCurve) = curve.label == SYSTEMATIC_TREND_LABEL
_kind(curve::ExtractedCurve) = _is_trend(curve) ? "systematic_trend" : "dataset"

# The variance of a linear functional `wᵀ R_T` of the curve from its propagated covariance alone.
function _propagated_variance(curve::ExtractedCurve, w::AbstractVector{<:Real})
    return max(dot(w, curve.R_T_covariance, w), 0.0)
end

# The variance of a linear functional `wᵀ R_T` of the curve: the propagated one, and for a curve
# that carries a jackknife the larger of that and the jackknife of the functional.
function _functional_variance(curve::ExtractedCurve, w::AbstractVector{<:Real})
    propagated = _propagated_variance(curve, w)
    Σ_J = curve.R_T_jackknife_covariance
    Σ_J === nothing && return propagated
    return max(propagated, dot(w, Σ_J, w))
end

"""
    range_mean(curve) -> Tuple{Float64,Float64}

The unweighted mean of the temperature ratio of a segmented curve over its mass range, with the
uncertainty propagated through the curve's covariance. For a systematic trend that carries a
jackknife, the term from the curve is the larger of the propagated one and the jackknife of the
same average over the leave-one-out refits.

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
rather than the sum of squares the independent-points form takes. For a systematic trend that
carries a jackknife, the term from the curve is the larger of the propagated one and the
jackknife of the same average over the leave-one-out refits. The yield term is unchanged.
The independent-points uncertainty, which is the published approximation and remains available
as `total_average(curve.R_T, yields)`, understates the propagated one; the method page
quantifies the difference for the shipped runs.

Throws as [`total_average`](@ref) does when the curve and the distribution share no mass number.
"""
function FissionFragmentsDomain.total_average(curve::ExtractedCurve, yields::MassYield)
    return _total_average(curve, yields, _functional_variance)
end

# The total average of a curve with the term from the curve given by `curve_variance(curve, w)`.
function _total_average(curve::ExtractedCurve, yields::MassYield, curve_variance)
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
    variance = curve_variance(curve, w)
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
[`yield_fraction`](@ref). For a systematic trend that carries a jackknife, the term from the
curve is the larger of the propagated one and the jackknife of the same average over the
leave-one-out refits; `uncertainty_propagated` is the uncertainty with the propagated term
alone, and equals `uncertainty` for every other curve.

    TotalAverage(curve, yields, heavy_masses)

The total average of `curve` over `yields`, `heavy_masses` being the heavy mass numbers of the
fragmentation range.
"""
struct TotalAverage
    value::Float64
    uncertainty::Float64
    uncertainty_independent_points::Float64
    yield_fraction::Float64
    uncertainty_propagated::Float64
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
        last(_total_average(curve, yields, _propagated_variance)),
    )
end
