# Method

## Temperature ratio from multiplicity data

Two premises carry the extraction. The prompt neutron multiplicity ratio of complementary
fragments follows their excitation energy ratio, and the fragments are excited highly enough for
their level densities to be of Fermi-gas form, ``E^* = a T^2``. Together,

```math
\frac{E_L^*}{E_H^*} = \frac{a_L T_L^2}{a_H T_H^2} \approx \frac{\nu_L}{\nu_H},
```

so with ``r_\nu = \nu_H/(\nu_L + \nu_H)`` and ``R_a = a_L/a_H``,

```math
R_T = \left[\frac{1 - r_\nu}{R_a\, r_\nu}\right]^{1/2}.
```

No fit and no prompt emission calculation enters. What remains to be chosen is the fragmentation
range, the prescription for the level density parameter, and the order in which the parameter
ratio is averaged over the isobaric charge distribution.

## Exact behaviour at the symmetric split

For a fissioning nucleus of even mass number the two fragments of the symmetric split are the same
nuclide. Three consequences are imposed by construction rather than fitted:

- the charge polarization vanishes, ``\Delta Z(A_0/2) = 0``, since there is nothing to polarize;
- the level density parameter ratio is exactly one, provided the retained charges are invariant
  under ``Z \to Z_0 - Z``, which the vanishing polarization secures;
- the multiplicity ratio is exactly one half, and therefore ``R_T(A_0/2) = 1``.

The last of these is not imposed on ``R_T``: it follows from pinning ``r_\nu`` and from
``R_a = 1``. A run records whether the identities hold, so that an input violating them is
reported rather than silently absorbed.

This also distinguishes the two averaging orders. [`RatioOfMeans`](@ref) satisfies
``R_a(A_0/2) = 1`` exactly. [`MeanOfRatios`](@ref) cannot: at the symmetric split the per-charge
ratios come in reciprocal pairs of equal weight, and the arithmetic mean of ``x`` and ``1/x``
exceeds one unless ``x = 1``. The former is the default for that reason.

## Parameterization by joined segments

The ratio extracted point by point is scattered, and for some data sets sparse, so it is
``r_\nu`` that is parameterized and ``R_T`` that follows by the exact relation above. The model is
continuous and piecewise-linear, written in the truncated-power basis

```math
f(x) = \beta_0 + \beta_1 (x - x_0) + \sum_k \gamma_k (x - \psi_k)_+,
```

in which continuity at every breakpoint holds identically and the fit at fixed breakpoints is one
weighted linear least-squares solve, yielding the parameter covariance. Breakpoints are restricted
to abscissae present in the data and searched exhaustively, which is tractable over a fragment
mass range and returns the global optimum. The number of segments is chosen by the Bayesian
information criterion, which prices both the extra slope and the extra breakpoint; the criterion
for every order examined is retained, so the choice can be audited.

Two constraints are enforced during the search because they are exact, not preferences: the fit
is pinned at the symmetric split, and it may not leave the interval ``(0, 1)``, outside which the
temperature ratio relation is undefined. Since a piecewise-linear function attains its extrema at
its pivots, the bound is tested exactly rather than sampled. Features that are strongly supported
but not exact — the minimum near ``A_H = 130``, the near-linear rise above the most probable
fragmentation — are left to the data, and can be imposed through `required_windows` when a sparse
data set does not determine them.

## References

- Eur. Phys. J. A **60**, 190 (2024) — the method and its conventions.
- T. von Egidy, D. Bucurescu, Phys. Rev. C **80**, 054310 (2009) — level density systematics.
- A. C. Wahl, At. Data Nucl. Data Tables **38**, 1 (1988) — charge polarization and dispersion.
- V. M. R. Muggeo, Stat. Med. **22**, 3055 (2003) — regression with unknown breakpoints.
- G. Schwarz, Ann. Stat. **6**, 461 (1978) — the information criterion.
