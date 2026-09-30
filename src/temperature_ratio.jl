# The temperature ratio of complementary fragments from a tabulated multiplicity ratio.
#
# The relation itself, `r_ν = E*_H/TXE` reduced over the isobaric charge distribution, its inverse
# and its slope, are those of FissionFragmentsDomain, which the consuming emission code applies in
# the forward direction. Only the tabulation over a ratio curve is added here.

"""
    temperature_ratio(averaging, model, domain, r_ν::RatioCurve) -> RatioCurve

The temperature ratio of complementary fully accelerated fragments at every heavy mass number of
the multiplicity ratio curve `r_ν = ν_H/(ν_L + ν_H)`, obtained by inverting the relation between
`R_T` and `E*_H/TXE` at each `A_H` of `domain` with the level density parameters of `model`, reduced
over the charge distribution as `averaging` prescribes:

- `ChargeResolved()` solves `r_ν = Σ_Z p(Z, A_H) / (1 + ρ_Z R_T²) / Σ_Z p(Z, A_H)`, with
  `ρ_Z = a_L/a_H` of each fragmentation, and `ChargeResolved(excitation)` weights each term by the
  mean total excitation of its fragmentation;
- `RatioOfMeans()` and `MeanOfRatios()` take the closed form `R_T = [(1 - r_ν)/(R_a r_ν)]^(1/2)`
  of Eur. Phys. J. A 60, 190 (2024), eq. (4), with the effective ratio `R_a = ⟨a_L⟩/⟨a_H⟩` or
  `⟨a_L/a_H⟩`.

The first is the exact inverse of a partition that applies `R_T` to every fragment pair with its
own parameters; see the method page of the documentation.

Uncertainties are carried from those of `r_ν` alone by the exact slope `∂R_T/∂r_ν` of the same
relation, `temperature_ratio_slope`; the level density systematics supply none. Where `r_ν` quotes
no uncertainty the result quotes none.

Mass numbers where the relation is undefined — outside the domain, no fragmentation with both
level density parameters, or `r_ν` outside `(0, 1)` — are omitted rather than extrapolated.
"""
function FissionFragmentsDomain.temperature_ratio(
    averaging::RatioAveraging,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    r_ν::RatioCurve,
)
    A_H = Int[]
    R_T = Float64[]
    σ = Union{Missing,Float64}[]

    for (index, mass) in enumerate(r_ν.A_H)
        mass in domain.heavy_masses || continue
        r = r_ν.ratio[index]
        0 < r < 1 || continue
        value = temperature_ratio(averaging, model, domain, mass, r)
        value === nothing && continue
        push!(A_H, mass)
        push!(R_T, value)
        push!(σ, _propagate(r_ν.σ[index], averaging, model, domain, mass, value))
    end

    return RatioCurve(A_H, R_T, σ, r_ν.label)
end

_propagate(::Missing, averaging, model, domain, A_H, R_T) = missing
function _propagate(σ_r::Real, averaging, model, domain, A_H, R_T)
    # Defined wherever the inverse is, which the caller has established.
    slope = temperature_ratio_slope(averaging, model, domain, A_H, R_T)
    return slope === nothing ? missing : abs(slope) * σ_r
end
