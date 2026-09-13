# Shared test inputs.
#
# The input data is held locally and is not shipped with the package, so a clone has none of it.
# Tests that exercise the real systematics on real nuclides are skipped when it is absent, and the
# rest of the suite — which uses synthetic inputs — runs regardless. `data/README.md` records what
# the files are and where they come from.

using Statistics: median

const DATA_DIRECTORY = joinpath(pkgdir(FissionTemperatureRatio), "data")
const MASS_EXCESS_FILE = joinpath(DATA_DIRECTORY, "reference", "mass_excess_ame2020.dat")
const SHELL_CORRECTION_FILE = joinpath(
    DATA_DIRECTORY, "reference", "shell_corrections_gilbert_cameron.dat"
)

const DATA_AVAILABLE = isfile(MASS_EXCESS_FILE) && isfile(SHELL_CORRECTION_FILE)

DATA_AVAILABLE || @warn "input data not present; the tests that need it are skipped" directory =
    DATA_DIRECTORY

const TEST_MASSES = DATA_AVAILABLE ? read_mass_excess_table(MASS_EXCESS_FILE) : nothing
const TEST_SHELLS =
    DATA_AVAILABLE ? read_shell_correction_table(SHELL_CORRECTION_FILE) : nothing
const BSFG_MODEL = DATA_AVAILABLE ? BackShiftedFermiGas(TEST_MASSES) : nothing
const GC_MODEL = DATA_AVAILABLE ? GilbertCameron(TEST_SHELLS) : nothing

const FLAT_CHARGES = ChargeDistribution(
    Dict{Int,Float64}(), Dict{Int,Float64}(), -0.5, 0.6, nothing
)

"A piecewise-linear reference shaped like the multiplicity ratio: a fall to a minimum, then a rise."
reference_ratio(A_H) = A_H ≤ 130 ? 0.5 - 0.03 * (A_H - 126) : 0.38 + 0.0125 * (A_H - 130)
