# Level density parameter systematics for nuclei occurring as fission fragments.

"""
    LevelDensityModel

Abstract supertype for the models that supply the level density parameter `a(A, Z)` of a
fission fragment.

Only energy-independent models are admissible here: the temperature ratio is extracted
from the Fermi-gas relation `E* = a T²` evaluated at a single, unknown excitation energy, so a
parameter that itself depends on `E*` would make the relation implicit.

A model carries the tabulated data it is evaluated from, so that
[`level_density_parameter`](@ref) needs nothing beyond the nuclide.
"""
abstract type LevelDensityModel end

"""
    BackShiftedFermiGas <: LevelDensityModel

Back-shifted Fermi gas systematic of von Egidy and Bucurescu, Phys. Rev. C **72**, 044311 (2005),
doi:10.1103/PhysRevC.72.044311, erratum Phys. Rev. C **73**, 049901 (2006),
doi:10.1103/PhysRevC.73.049901, and Phys. Rev. C **80**, 054310 (2009),
doi:10.1103/PhysRevC.80.054310,

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

For the majority of nuclei occurring as fission fragments this systematic reproduces the
superfluid-model level density parameter at the excitation energies fragments actually attain,
which is why it is the default here.
"""
struct BackShiftedFermiGas <: LevelDensityModel
    masses::MassExcessTable
end

"""
    LiquidDropCoefficients

Coefficients of the liquid-drop binding energy entering the shell correction of
[`BackShiftedFermiGas`](@ref): volume, surface, Coulomb, and the two coefficients of the
mass-dependent symmetry term `a_sym = A (c₁ - c₂ A^(-1/3))`. Values in MeV, from the liquid-drop
mass formula of J. M. Pearson, Hyperfine Interact. **132**, 59 (2001),
doi:10.1023/A:1011973100463, as adopted in Phys. Rev. C **72**, 044311 (2005).
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
P_d = [Δ(A+2, Z+1) - 2 Δ(A, Z) + Δ(A-2, Z-1)] / 4.
```

Returns `missing` when any of the three mass excesses required for `P_d` is absent from the
table, which is the normal situation at the edges of the known region.
"""
function shell_correction(A::Integer, Z::Integer, masses::MassExcessTable)
    Δ = mass_excess(masses, A, Z)
    Δ₊ = mass_excess(masses, A + 2, Z + 1)
    Δ₋ = mass_excess(masses, A - 2, Z - 1)
    (ismissing(Δ) || ismissing(Δ₊) || ismissing(Δ₋)) && return missing

    W_exp = Z * masses.Δᵖ + (A - Z) * masses.Δⁿ - Δ
    η = (A - 2 * Z) / A
    a_sym = A * (LIQUID_DROP.symmetry_constant - LIQUID_DROP.symmetry_slope * A^(-1 / 3))
    W_LDM =
        LIQUID_DROP.volume * A - LIQUID_DROP.surface * A^(2 / 3) -
        LIQUID_DROP.coulomb * Z^2 * A^(-1 / 3) - a_sym * η^2
    P_d = (Δ₊ - 2 * Δ + Δ₋) / 4

    return (W_LDM - W_exp) + P_d
end

"""
    level_density_parameter(model, A, Z) -> Union{Float64,Missing}

Level density parameter `a` in MeV⁻¹ of the nuclide `(A, Z)`.

Returns `missing` when the model cannot be evaluated — the tabulated data it needs is
absent for this nuclide, or the expression yields a non-positive value, which is not physically
meaningful.

# Example

```julia
masses = read_mass_excess_table(joinpath(datadir(), "reference", "mass_excess_ame2020.dat"))
a = level_density_parameter(BackShiftedFermiGas(masses), 132, 50)
```
"""
function level_density_parameter(model::BackShiftedFermiGas, A::Integer, Z::Integer)
    δW = shell_correction(A, Z, model.masses)
    ismissing(δW) && return missing
    a = (BSFG.p₁ + BSFG.p₂ * δW) * A^BSFG.p₃
    return a > 0 ? a : missing
end

"""
    GilbertCameron <: LevelDensityModel

Level density systematic of Gilbert and Cameron, Can. J. Phys. **43**, 1446 (1965),
doi:10.1139/p65-139, Eq. (20),

```
a = A [c₁ (S_Z + S_N) + c₂],
```

with `S_Z` and `S_N` the shell corrections of that paper's Table III, tabulated against proton and
neutron number.

Eq. (20) is the correlation the authors fit to **undeformed** nuclei. They give a second line,
parallel to it, for deformed nuclei — Eq. (21), the same slope with an offset of `0.120` in place
of `0.142`, some fifteen per cent lower. Only Eq. (20) is implemented. That is deliberate: the
data behind the deformed branch is thin, and the resulting uncertainty in where the deformation
boundary falls is part of why this systematic was later superseded by the back-shifted Fermi gas
of [`BackShiftedFermiGas`](@ref), which needs no such division. Applying Eq. (20) throughout is
therefore the same choice the published results of this method were obtained under, and it is
retained for that reason rather than for want of the alternative.

It is provided for assessing how much the extracted temperature ratio depends on the level density
model, not as an equal alternative: for nuclei occurring as fission fragments it returns
parameters well above those of the superfluid model and of the back-shifted Fermi gas, so the
spread between the two models bounds a systematic uncertainty the propagated experimental
uncertainties do not cover.
"""
struct GilbertCameron <: LevelDensityModel
    shells::ShellCorrectionTable
end

"""
    GilbertCameronCoefficients

The two coefficients of `a/A = c₁ (S_Z + S_N) + c₂` of Gilbert and Cameron, Can. J. Phys. **43**,
1446 (1965), Eq. (20), the correlation fitted to undeformed nuclei. Their Eq. (21) keeps `c₁` and
replaces `c₂` by `0.120` for deformed nuclei; see [`GilbertCameron`](@ref).
"""
struct GilbertCameronCoefficients
    c₁::Float64
    c₂::Float64
end

const GILBERT_CAMERON = GilbertCameronCoefficients(9.17e-3, 1.42e-1)

function level_density_parameter(model::GilbertCameron, A::Integer, Z::Integer)
    N = A - Z
    S_Z = get(model.shells.S_Z, Int(Z), missing)
    S_N = get(model.shells.S_N, Int(N), missing)
    (ismissing(S_Z) || ismissing(S_N)) && return missing
    a = A * (GILBERT_CAMERON.c₁ * (S_Z + S_N) + GILBERT_CAMERON.c₂)
    return a > 0 ? a : missing
end
