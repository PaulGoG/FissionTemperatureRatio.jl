# Level density parameter systematics for nuclei occurring as fission fragments.

"""
    LevelDensityPrescription

Abstract supertype for the prescriptions that supply the level density parameter `a(A, Z)` of a
fission fragment.

Only energy-independent prescriptions are admissible here: the temperature ratio is extracted
from the Fermi-gas relation `E* = a T²` evaluated at a single, unknown excitation energy, so a
parameter that itself depends on `E*` would make the relation implicit.

A prescription carries the tabulated data it is evaluated from, so that
[`level_density_parameter`](@ref) needs nothing beyond the nuclide.
"""
abstract type LevelDensityPrescription end

"""
    BackShiftedFermiGas <: LevelDensityPrescription

Back-shifted Fermi gas systematic of von Egidy and Bucurescu, Phys. Rev. C **72**, 044311 (2005),
erratum Phys. Rev. C **73**, 049901 (2006), and Phys. Rev. C **80**, 054310 (2009),

```
a = (p₁ + p₂ δW) A^p₃,   δW = δW₀ + P_d,
```

with `δW₀` the shell correction and `P_d` the deuteron pairing energy.

The shell correction is `S(Z,N) = M_exp - M_LD` of Phys. Rev. C **72**, 044311 (2005), Eq. (7), a
difference of masses; since `M = Z M_p + N M_n - E_b`, it is equivalently the difference of the
liquid-drop and experimental binding energies, which is how it is computed here. The pairing term
is `0.5 P'a` with `P'a` of Phys. Rev. C **80**, 054310 (2009), Eq. (12); the 2009 paper absorbs
the sign alternation of the earlier definition into `P'a`, so no structure-dependent case
distinction is needed.

Only `a` is taken from this systematics. The papers fit it jointly with the back-shift
`E1 = -0.381 + 0.5 P'a` of Phys. Rev. C **80**, 054310 (2009), Eq. (20), for the level density
`ρ(U) ∝ exp(2√(a(U - E1)))`, whose closure is `U - E1 = a T²`. The extraction implemented here
rests instead on the un-shifted `E* = a T²`, which is the premise the method is published under.
The shift does not cancel in the ratio, since `E1` depends on the fragment through `P'a` and
`E1_L ≠ E1_H`. It is of order ±1 MeV against fragment excitations of 10-20 MeV, and its effect on
`R_T` has not been quantified.

For the majority of nuclei occurring as fission
fragments this systematic reproduces the superfluid-model level density parameter at the
excitation energies fragments actually attain, which is why it is the default here.
"""
struct BackShiftedFermiGas <: LevelDensityPrescription
    masses::MassExcessTable
end

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
`a = (p₁ + p₂ δW) A^p₃` of von Egidy and Bucurescu, Phys. Rev. C **80**, 054310 (2009), Eq. (19):
`p₁ = 0.199(7)`, `p₂ = 0.0096(4)`, `p₃ = 0.869(7)`.
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
    level_density_parameter(prescription, A, Z) -> Union{Float64,Missing}

Level density parameter `a` in MeV⁻¹ of the nuclide `(A, Z)`.

Returns `missing` when the prescription cannot be evaluated — the tabulated data it needs is
absent for this nuclide, or the expression yields a non-positive value, which is not physically
meaningful.

# Example

```julia
masses = read_mass_excess(joinpath(datadir(), "mass_excess", "AME2020.ANA"))
a = level_density_parameter(BackShiftedFermiGas(masses), 132, 50)
```
"""
function level_density_parameter(prescription::BackShiftedFermiGas, A::Integer, Z::Integer)
    δW = shell_correction(A, Z, prescription.masses)
    ismissing(δW) && return missing
    a = (BSFG.p₁ + BSFG.p₂ * δW) * A^BSFG.p₃
    return a > 0 ? a : missing
end

"""
    GilbertCameron <: LevelDensityPrescription

Level density systematic of Gilbert and Cameron, Can. J. Phys. **43**, 1446 (1965), for spherical
nuclei,

```
a = A [c₁ (S_Z + S_N) + c₂],
```

with `S_Z` and `S_N` the tabulated shell corrections.

It is provided for assessing how much the extracted temperature ratio depends on the level density
prescription, not as an equal alternative: for nuclei occurring as fission fragments it returns
parameters well above those of the superfluid model and of the back-shifted Fermi gas, so the
spread between the two prescriptions bounds a systematic uncertainty the propagated experimental
uncertainties do not cover.
"""
struct GilbertCameron <: LevelDensityPrescription
    shells::ShellCorrectionTable
end

"""
    GilbertCameronCoefficients

The two coefficients of `a/A = c₁ (S_Z + S_N) + c₂` of Gilbert and Cameron, Can. J. Phys. **43**,
1446 (1965).
"""
struct GilbertCameronCoefficients
    c₁::Float64
    c₂::Float64
end

const GILBERT_CAMERON = GilbertCameronCoefficients(9.17e-3, 1.42e-1)

function level_density_parameter(prescription::GilbertCameron, A::Integer, Z::Integer)
    N = A - Z
    S_Z = get(prescription.shells.S_Z, Int(Z), missing)
    S_N = get(prescription.shells.S_N, Int(N), missing)
    (ismissing(S_Z) || ismissing(S_N)) && return missing
    a = A * (GILBERT_CAMERON.c₁ * (S_Z + S_N) + GILBERT_CAMERON.c₂)
    return a > 0 ? a : missing
end
