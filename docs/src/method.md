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
range, the level density model, and the order in which the parameter ratio is averaged over the
isobaric charge distribution.

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

## Description by joined segments

The ratio extracted point by point is scattered, and for some datasets sparse, so it is
``r_\nu`` that is described by joined segments and ``R_T`` that follows by the exact relation above. The model is
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
is pinned to one half at the symmetric split — for a dataset whose complete pairs begin above
it the identity has no abscissa to act on, and that curve is fitted unpinned — and it may not leave the interval ``(0, 1)``, outside which the
temperature ratio relation is undefined. Since a piecewise-linear function attains its extrema at
its pivots, the bound is tested exactly rather than sampled.

## One curve per dataset, and one systematic trend

The datasets of a fissioning nucleus can differ well beyond their quoted uncertainties, and where
they do, the temperature ratio can only be determined separately for each of them. A run therefore
produces one segmented curve per dataset, fitted to that dataset alone.

Alongside them it produces a systematic-trend curve, fitted through the whole body of data with
the minimum at the heavy magic fragment *placed* rather than fitted — `required_windows` in the
configuration, defaulting to ``A_H \in [128, 132]``. Its rise above the most probable
fragmentation therefore falls between those of the individual sets. This is the curve to use where
a dataset is too sparse or too scattered to resolve the shape on its own, and where no energy
partition from a scission model is available.

These are alternatives, not an ensemble to be averaged. A prompt emission code takes one of them
as input; which one describes reality is settled downstream, by comparing the prompt neutron
multiplicity distributions and the fragment yields that code produces against experimental data.
The package's job is to supply the candidates, each traceable to the measurement it came from.

## The total average

Prompt emission codes that take a single temperature ratio for all fragmentations, rather than a
function of mass number, need one number. It is obtained by averaging over a fragment mass yield
distribution,

```math
\langle R_T \rangle = \frac{\sum_{A_H} Y(A_H)\, R_T(A_H)}{\sum_{A_H} Y(A_H)},
```

taken over the heavy branch, with ``Y`` the pre-neutron mass yield — the temperature ratio is a
function of the primary heavy-fragment mass number, so a post-neutron distribution would weight
each ratio by the yield of a different fragmentation. The normalization of ``Y`` cancels.

This is the quantity the literature tabulates, and it is not the mean over the fragment mass
range. That mean weights every mass number equally, so the far-asymmetric tail, where the yield is
smaller by orders of magnitude, counts as much as the peak. Both are reported, under names that
distinguish them.

Uncertainties propagate from the ratio and from the yield,

```math
\sigma^2 = \sum_{A_H} \left[ \left(\frac{Y}{\sum Y}\right)^2 \sigma_{R_T}^2
         + \left(\frac{R_T - \langle R_T \rangle}{\sum Y}\right)^2 \sigma_Y^2 \right],
```

which is why a multiplicity dataset quoting no uncertainties still yields an uncertain average:
the yield distribution supplies it. Correlations between mass numbers are neglected in both
inputs, the sources not reporting them.

## Combining datasets

Several measurements of one fissioning system are not merged into one. Each is fitted on its own,
and a further curve is fitted to their combination.

That combination is not a concatenation. At each mass number the available values are combined
with inverse-variance weights carrying an additional between-dataset variance ``\tau^2``,
estimated from their dispersion after DerSimonian and Laird,

```math
w = \frac{1}{\sigma^2 + \tau^2}, \qquad
\bar{r} = \frac{\sum w\, r}{\sum w}, \qquad
\sigma_{\bar{r}} = \left(\sum w\right)^{-1/2}.
```

The reason is empirical. The datasets of one system disagree by ten to twenty times their quoted
uncertainties, so ``\tau^2`` dominates ``\sigma^2``, the weights become nearly equal, and the
combination stops being decided by whichever author quoted the smallest errors. It also makes the
fitted chi-squared of the combined curve a statement about the fit rather than about the
disagreement: for 252-Cf it falls from about 16 to about 2.

The same disagreement is why no dataset is rejected for being far from the others. In units of
the quoted uncertainties none of them agrees with any other, so such a criterion rejects whatever
it is tuned to reject. What can be said about a dataset without reference to the rest — how many
usable fragment pairs it has, whether its ratio stays inside ``(0,1)``, whether it satisfies the
identity at the symmetric split, whether complementary multiplicities sum to the total — is
reported for every dataset, and exclusions are named explicitly rather than inferred.

## References

- Eur. Phys. J. A **60**, 190 (2024) — the method and its conventions.
- T. von Egidy, D. Bucurescu, Phys. Rev. C **80**, 054310 (2009) — level density systematics.
- A. C. Wahl, At. Data Nucl. Data Tables **38**, 1 (1988) — charge polarization and dispersion.
- V. M. R. Muggeo, Stat. Med. **22**, 3055 (2003) — regression with unknown breakpoints.
- G. Schwarz, Ann. Stat. **6**, 461 (1978) — the information criterion.
- R. DerSimonian, N. Laird, Control. Clin. Trials **7**, 177 (1986) — the between-set variance.
