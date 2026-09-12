# Level density parameter systematics for nuclei occurring as fission fragments.

"""
    LevelDensityPrescription

Abstract supertype for the prescriptions that supply the level density parameter `a(A, Z)` of a
fission fragment.

Only energy-independent prescriptions are admissible here: the temperature ratio is extracted
from the Fermi-gas relation `E* = a T²` evaluated at a single, unknown excitation energy, so a
parameter that itself depends on `E*` would make the relation implicit.
"""
abstract type LevelDensityPrescription end

"""
    BackShiftedFermiGas <: LevelDensityPrescription

Back-shifted Fermi gas systematic of von Egidy and Bucurescu, Phys. Rev. C **80**, 054310 (2009),

```
a = (p₁ + p₂ δW) A^p₃,   δW = δW₀ + P_d,
```

with `δW₀` the shell correction, the difference between the liquid-drop and experimental binding
energies, and `P_d` the deuteron pairing energy. For the majority of nuclei occurring as fission
fragments this systematic reproduces the superfluid-model level density parameter at the
excitation energies fragments actually attain, which is why it is the default here.
"""
struct BackShiftedFermiGas <: LevelDensityPrescription end

"""
    LiquidDropCoefficients

Coefficients of the liquid-drop binding energy entering the shell correction of
[`BackShiftedFermiGas`](@ref): volume, surface, Coulomb, and the two coefficients of the
mass-dependent symmetry term `a_sym = A (c₁ - c₂ A^(-1/3))`. Values in MeV, after Myers and
Swiatecki.
"""
struct LiquidDropCoefficients
    volume::Float64
    surface::Float64
    coulomb::Float64
    symmetry_constant::Float64
    symmetry_slope::Float64
end

const LIQUID_DROP = LiquidDropCoefficients(15.65, 17.63, 0.864 / 1.233, 27.72, 25.6)

"""
    BSFGCoefficients

The three fitted coefficients of the back-shifted Fermi gas expression
`a = (p₁ + p₂ δW) A^p₃` of von Egidy and Bucurescu, Phys. Rev. C **80**, 054310 (2009).
"""
struct BSFGCoefficients
    p₁::Float64
    p₂::Float64
    p₃::Float64
end

const BSFG = BSFGCoefficients(1.99e-1, 9.6e-3, 8.69e-1)

"""
    shell_correction(A, Z, masses) -> Union{Float64,Missing}

Shell correction `δW = δW₀ + P_d` in MeV for the nuclide `(A, Z)`.

`δW₀` is the difference between the liquid-drop binding energy and the experimental one obtained
from the mass excess table `masses`; `P_d` is the deuteron pairing energy, the second difference
of the mass excess along the line of constant neutron excess,

```
P_d = [D(A+2, Z+1) - 2 D(A, Z) + D(A-2, Z-1)] / 4.
```

Returns `missing` when any of the three mass excesses required for `P_d` is absent from the
table, which is the normal situation at the edges of the known region.
"""
function shell_correction(A::Integer, Z::Integer, masses::MassExcessTable)
    D = mass_excess(masses, A, Z)
    D₊ = mass_excess(masses, A + 2, Z + 1)
    D₋ = mass_excess(masses, A - 2, Z - 1)
    (ismissing(D) || ismissing(D₊) || ismissing(D₋)) && return missing

    W_exp = Z * masses.Dᵖ + (A - Z) * masses.Dⁿ - D
    η = (A - 2 * Z) / A
    a_sym = A * (LIQUID_DROP.symmetry_constant - LIQUID_DROP.symmetry_slope * A^(-1 / 3))
    W_LDM =
        LIQUID_DROP.volume * A - LIQUID_DROP.surface * A^(2 / 3) -
        LIQUID_DROP.coulomb * Z^2 * A^(-1 / 3) - a_sym * η^2
    P_d = (D₊ - 2 * D + D₋) / 4

    return (W_LDM - W_exp) + P_d
end

"""
    level_density_parameter(prescription, A, Z, masses) -> Union{Float64,Missing}

Level density parameter `a` in MeV⁻¹ of the nuclide `(A, Z)`.

Returns `missing` when the prescription cannot be evaluated — the required mass excesses are
absent, or the expression yields a non-positive value, which happens for a few nuclides with
large negative shell corrections and is not physically meaningful.

# Example

```jldoctest
julia> masses = read_mass_excess(joinpath(datadir(), "mass_excess", "AME2020.ANA"));

julia> a = level_density_parameter(BackShiftedFermiGas(), 132, 50, masses);

julia> 10 < a < 20
true
```
"""
function level_density_parameter(
    ::BackShiftedFermiGas, A::Integer, Z::Integer, masses::MassExcessTable
)
    δW = shell_correction(A, Z, masses)
    ismissing(δW) && return missing
    a = (BSFG.p₁ + BSFG.p₂ * δW) * A^BSFG.p₃
    return a > 0 ? a : missing
end
