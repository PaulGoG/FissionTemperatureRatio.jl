# Shared test inputs. The mass excesses are the table the package ships, so that the level density
# tests exercise the real systematic on real nuclides rather than a surrogate that would only
# re-derive the formula under test.

using Statistics: median

const TEST_MASSES = read_mass_excess(
    joinpath(pkgdir(FissionTemperatureRatio), "data", "mass_excess", "AME2020.ANA")
)

const FLAT_CHARGES = ChargeDistributionData(
    Dict{Int,Float64}(), Dict{Int,Float64}(), -0.5, 0.6, nothing
)

"A piecewise-linear reference shaped like the multiplicity ratio: a fall to a minimum, then a rise."
reference_ratio(A_H) = A_H ≤ 130 ? 0.5 - 0.03 * (A_H - 126) : 0.38 + 0.0125 * (A_H - 130)
