# Level density parameter ratio and the temperature ratio of complementary fragments.

"""
    RatioAveraging

Abstract supertype for the two orders in which the level density parameter ratio of complementary
fragments can be reduced over the isobaric charge distribution.

The two differ by Jensen's inequality and coincide only for a degenerate charge distribution.
Reported differences over the fragment mass range are small, and largest around the heavy magic
fragment, which is also where the multiplicity ratio has its minimum — so the choice is recorded
with every result rather than left implicit.
"""
abstract type RatioAveraging end

"""
    RatioOfMeans <: RatioAveraging

Average the level density parameters and then form their ratio, `R_a = ⟨a_L⟩ / ⟨a_H⟩`.

This asserts the Fermi-gas relation for the mass-resolved fragment, `⟨E*⟩ = ⟨a⟩ T²`, with a
temperature that is a function of mass number alone, which is the natural order when the measured
input, `ν(A)`, is itself resolved by mass number only.

It is the default because it alone respects the exact behaviour at the symmetric split. There the
two fragments are the same nuclide and the charges retained are invariant under `Z -> Z₀ - Z`, so
the two averages run over the same set and their ratio is exactly one, giving `R_T(A₀/2) = 1` as
an outcome rather than an imposition. See [`MeanOfRatios`](@ref) for why the other order cannot
do this.
"""
struct RatioOfMeans <: RatioAveraging end

"""
    MeanOfRatios <: RatioAveraging

Form the ratio at each charge and then average it, `R_a = ⟨a_L / a_H⟩`.

This is the order obtained by writing the excitation energy partition for each fragmentation
`(A_H, Z_H)` separately and averaging afterwards. To first order in the dispersion of the ratio
over the charge distribution it is the effective ratio that inverts a charge-averaged
multiplicity ratio, which is an argument in its favour away from symmetry.

At the symmetric split, however, it is provably wrong. The per-charge ratios there come in
reciprocal pairs of equal weight, and the arithmetic mean of `x` and `1/x` exceeds one unless
`x = 1`, so this order returns `R_a > 1` where the two fragments are the same nuclide and the
ratio must be exactly one. The departure is small — of order a part in a thousand — but it is
systematic and one-sided, and it propagates into `R_T` at the one mass number where the answer is
known in advance. Offered for comparison; not the default.
"""
struct MeanOfRatios <: RatioAveraging end

"""
    level_density_ratio(averaging, prescription, A₀, Z₀, domain) -> Dict{Int,Float64}

The level density parameter ratio of complementary fragments,

```
R_a(A_H) = a_L / a_H,
```

averaged over the isobaric charge distribution in the order fixed by `averaging`, for every heavy
mass number of the domain's range.

Mass numbers at which neither fragment of any charge pair has a defined level density parameter
are absent from the result.
"""
function level_density_ratio(
    averaging::RatioAveraging,
    prescription::LevelDensityPrescription,
    A₀::Integer,
    Z₀::Integer,
    domain::FragmentationDomain,
)
    ratio = Dict{Int,Float64}()
    for A_H in domain.A_H_range
        A_L = A₀ - A_H
        value = _charge_averaged_ratio(averaging, prescription, A_H, A_L, Z₀, domain)
        ismissing(value) || value ≤ 0 || (ratio[A_H] = value)
    end
    return ratio
end

function _charge_averaged_ratio(
    ::MeanOfRatios,
    prescription::LevelDensityPrescription,
    A_H::Integer,
    A_L::Integer,
    Z₀::Integer,
    domain::FragmentationDomain,
)
    ratios = Dict{Int,Union{Float64,Missing}}()
    for Z_H in charges(domain, A_H)
        a_H = level_density_parameter(prescription, A_H, Z_H)
        a_L = level_density_parameter(prescription, A_L, Z₀ - Z_H)
        ratios[Z_H] = (ismissing(a_H) || ismissing(a_L)) ? missing : a_L / a_H
    end
    return average_over_charge(ratios, domain, A_H)
end

function _charge_averaged_ratio(
    ::RatioOfMeans,
    prescription::LevelDensityPrescription,
    A_H::Integer,
    A_L::Integer,
    Z₀::Integer,
    domain::FragmentationDomain,
)
    heavy = Dict{Int,Union{Float64,Missing}}()
    light = Dict{Int,Union{Float64,Missing}}()
    for Z_H in charges(domain, A_H)
        a_H = level_density_parameter(prescription, A_H, Z_H)
        a_L = level_density_parameter(prescription, A_L, Z₀ - Z_H)
        # Both parameters of a pair must exist, so that the two averages run over the same charges.
        if ismissing(a_H) || ismissing(a_L)
            heavy[Z_H] = missing
            light[Z_H] = missing
        else
            heavy[Z_H] = a_H
            light[Z_H] = a_L
        end
    end
    ā_H = average_over_charge(heavy, domain, A_H)
    ā_L = average_over_charge(light, domain, A_H)
    (ismissing(ā_H) || ismissing(ā_L) || ā_H ≤ 0) && return missing
    return ā_L / ā_H
end

"""
    temperature_ratio(r_ν, R_a) -> RatioCurve

The temperature ratio of complementary fully accelerated fragments,

```
R_T(A_H) = T_L / T_H = [(1 - r_ν) / (R_a r_ν)]^(1/2),
```

obtained from the prompt neutron multiplicity ratio `r_ν = ν_H/(ν_L + ν_H)` and the level density
parameter ratio `R_a = a_L/a_H`.

The relation follows from the Fermi-gas form of the fragment excitation energies together with
the identification of their ratio with the multiplicity ratio, and it is exact given those two
premises — no fit and no prompt emission calculation enters.

Uncertainties are propagated from those of `r_ν` alone,

```
σ_RT = σ_r / (2 R_T R_a r_ν²),
```

since the level density parameter systematic supplies no uncertainty. The uncertainty of `R_T` is
therefore a lower bound: the spread between level density prescriptions is the larger effect, and
is assessed by repeating the extraction with [`GilbertCameron`](@ref) in place of
[`BackShiftedFermiGas`](@ref).

Mass numbers absent from `R_a` are omitted.
"""
function temperature_ratio(r_ν::RatioCurve, R_a::AbstractDict{Int,Float64})
    A_H = Int[]
    value = Float64[]
    σ = Float64[]

    for (index, mass) in enumerate(r_ν.A_H)
        haskey(R_a, mass) || continue
        ratio = R_a[mass]
        r = r_ν.value[index]
        (r > 0 && r < 1 && ratio > 0) || continue

        R_T = sqrt((1 - r) / (ratio * r))
        σ_R_T = r_ν.σ[index] / (2 * R_T * ratio * r^2)

        push!(A_H, mass)
        push!(value, R_T)
        push!(σ, σ_R_T)
    end

    return RatioCurve(A_H, value, σ, r_ν.label)
end

"""
    weighted_mean(curve) -> Tuple{Float64,Float64}

Inverse-variance weighted mean of a ratio curve and the uncertainty of that mean.

Points without an uncertainty carry no information about the weighting, so the mean falls back to
the unweighted one when no point in the curve has a positive uncertainty. The uncertainty of the
mean is `(Σ w)^(-1/2)` for the weighted case and the standard error of the mean otherwise; it is
not the quadrature sum of the input uncertainties, which would grow with the number of points.
"""
function weighted_mean(curve::RatioCurve)
    isempty(curve) && throw(ArgumentError("cannot average an empty ratio curve"))

    usable = curve.σ .> 0
    if !any(usable)
        n = length(curve)
        μ = sum(curve.value) / n
        n == 1 && return (μ, 0.0)
        return (μ, sqrt(sum(abs2, curve.value .- μ) / (n * (n - 1))))
    end

    w = 1 ./ curve.σ[usable] .^ 2
    μ = sum(w .* curve.value[usable]) / sum(w)
    return (μ, 1 / sqrt(sum(w)))
end
