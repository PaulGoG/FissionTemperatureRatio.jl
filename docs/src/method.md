# Method

## Temperature ratio from multiplicity data

The method and its conventions are those of [Tudora2024](@cite). Two premises carry the
extraction. The prompt neutron multiplicity ratio of complementary fragments follows their
excitation energy ratio, eq. (1) of that paper, and the fragments are excited highly enough for
their level densities to be of Fermi-gas form, ``E^* = a T^2``. For one fragmentation
``(A_H, Z)`` with ``\rho_Z = a_L/a_H``, the heavy fragment then carries

```math
\frac{E_H^*}{\mathrm{TXE}} = \frac{1}{1 + \rho_Z R_T^2}
```

of the total excitation energy. The measured multiplicity ratio ``r_\nu = \nu_H/(\nu_L + \nu_H)``
is resolved by mass alone: each ``\nu(A)`` is a yield-weighted mean over the charges and kinetic
energies of that mass. Carried to such means, the premise identifies ``r_\nu(A_H)`` with the
excitation-weighted mean of the relation over the isobaric charge distribution ``p(Z, A_H)``,

```math
r_\nu(A_H) = \frac{\sum_Z p(Z, A_H)\, \langle \mathrm{TXE} \rangle_Z \,/\, (1 + \rho_Z R_T^2)}
                  {\sum_Z p(Z, A_H)\, \langle \mathrm{TXE} \rangle_Z},
\qquad
\langle \mathrm{TXE} \rangle_Z = Q(A_H, Z) + E^*_{\mathrm{CN}} - \langle \mathrm{TKE} \rangle(A_H),
```

and ``R_T(A_H)`` is its root, which is unique: the right-hand side falls strictly from one to zero
as ``R_T`` grows. This is `ratio_averaging = "charge_resolved"`, the default. The mean total
kinetic energy comes from a measured pre-neutron ``\langle \mathrm{TKE} \rangle(A)`` named in the
configuration; without one, every fragmentation is weighted by ``p(Z, A_H)`` alone, and the run
says so in its log and in its manifest. A heavy mass beyond the measured span takes the value at
the nearest measured mass and is listed as extrapolated in the run metadata.

For a single effective ratio ``R_a`` in place of the per-charge ``\rho_Z`` the relation inverts in
closed form, eq. (4) of [Tudora2024](@cite),

```math
R_T = \left[\frac{1 - r_\nu}{R_a\, r_\nu}\right]^{1/2},
```

with ``R_a = \langle a_L \rangle / \langle a_H \rangle`` (`"ratio_of_means"`, the setting of the
published extraction) or ``R_a = \langle a_L / a_H \rangle`` (`"mean_of_ratios"`); both averages
run over ``p(Z, A_H)``. The three coincide for a single charge per mass and differ at second order
in the spread of ``\rho_Z`` over the charge window. The paper adopts the ratio of means and shows
it beside the mean of ratios in its Fig. 2; it remains selectable for that reason.

No fit and no prompt emission calculation enters. The fragmentation domain — heavy masses, the
charges retained about ``Z_p(A)``, their probabilities — the mass table, the level density
parameters and the relation between ``R_T`` and ``E_H^*/\mathrm{TXE}`` are those of
FissionFragmentsDomain.jl, which the codes that consume ``R_T(A_H)`` partition on as well, so the
two sides do not hold diverging copies of one definition. The charge distribution is Wahl's
``Z_p`` model with the parameters of [Wahl1988](@cite) for the four reactions evaluated there and
the systematics of [Wahl2002](@cite) otherwise. The level density parameter is taken from the
back-shifted Fermi-gas systematics of [vonEgidy2005, vonEgidy2009](@cite) by default, or from the
formula of [GilbertCameron1965](@cite) with the shell corrections of its Table III, whose
eq. (21), ``a/A = 0.00917\,S + 0.120``, applies to the nuclei it counts as deformed and eq. (20),
``0.00917\,S + 0.142``, to the rest; about a third of the yield-weighted ``^{252}``Cf fragments
lie in the first deformed region. The shell corrections ``S`` are referred to a spherical
reference, and eq. (21) is what corrects them for the deformation of those nuclei
[BrancazioCameron1969](@cite), so both formulas apply by default. `deformed_branch = false`
applies eq. (20) throughout, the spherical-only setting of Table 2 of [Tudora2024](@cite).

## Why the charge-resolved inversion

1. **It is the inverse of what the consumer does.** The Point-by-Point and sequential emission
   treatments apply a temperature ratio to every fragment pair ``(A, Z, \mathrm{TKE})`` with its
   own ``a_L/a_H``, eq. (5) of [Tudora2022](@cite), and so do the Monte Carlo codes CGMF
   [Talou2021](@cite) and FIFRELIN [Piau2023](@cite). The charge-resolved inversion closes that
   round trip to rounding. An effective ratio does not: on the ``^{252}``Cf(sf) domain with the
   back-shifted Fermi gas, the ratio of means misses it by ``2.9 \times 10^{-3}`` in ``R_T`` along
   the systematic trend and the mean of ratios by ``1.4 \times 10^{-2}``, both at the doubly magic
   heavy fragment, ``A_H = 132``; over ``0.8 \le R_T \le 1.45`` and the four shipped systems the
   misses reach ``4.2 \times 10^{-3}`` and ``1.3 \times 10^{-2}``.
2. **It returns ``R_T = 1`` from ``r_\nu = 1/2`` at the symmetric split** exactly, provided the
   charge set there is its own mirror under ``Z \to Z_0 - Z``. The terms of the sum then pair as
   ``\rho`` and ``1/\rho`` with equal weight, and ``1/(1 + \rho) + 1/(1 + 1/\rho) = 1``; the
   excitation weights pair the same way, ``Z`` and ``Z_0 - Z`` being one fragmentation.
3. **The excitation weight is the premise itself.** ``\nu_L/\nu_H = E_L^*/E_H^*`` holds, if at
   all, fragmentation by fragmentation; a mass-resolved multiplicity weights each fragmentation by
   the excitation it shares out. At the ``Z = 50`` shell the ``Q``-value and ``a_L/a_H`` change
   together across the charge window, so the weight matters there: omitting it moves ``R_T`` at
   ``A_H = 130`` by ``7 \times 10^{-3}`` to ``1.1 \times 10^{-2}`` along the systematic trends of
   the four shipped systems, with the ``\langle \mathrm{TKE} \rangle(A)`` of Göök et al.
   (``^{252}``Cf), Al-Adili et al. (``^{235}``U), Wagemans et al. (``^{239}``Pu) and Geltenbort
   et al. (``^{233}``U). ``\mathrm{TKE}`` is a property of the split, so where a dataset gives both
   ``A`` and ``A_0 - A`` the pair takes their mean. Double-energy measurements differ in their
   pulse-height-defect calibration by several MeV; a run records the offset of the input's
   yield-weighted mean from the energy standard of the system, the recommendation of Gönnenwein as
   tabulated in [Bertsch2015](@cite).
4. **FIFRELIN fixes ``R_T(A_{\mathrm{CN}}/2) = 1``** [Piau2023](@cite); the charge-resolved
   inversion reaches the same value from the data rather than by construction.
5. **Uncertainties propagate with the exact slope** ``\partial R_T/\partial r_\nu`` of the
   inversion, the reciprocal of
   ``-\sum_Z w_Z\, 2 \rho_Z R_T/(1 + \rho_Z R_T^2)^2`` with the normalized weights ``w_Z``,
   in place of the closed-form ``-1/(2 R_T R_a r_\nu^2)`` of an effective ratio.

## Exact behaviour at the symmetric split

For a fissioning nucleus of even mass number the two fragments of the symmetric split are the same
nuclide. The charge polarization vanishes there, ``\Delta Z(A_0/2) = 0``, as a property of the
charge model rather than an imposition: Wahl's 1988 fits pass through point X of their Fig. 2 at
the symmetric split, and the steep branch of the 2002 systematics crosses zero there linearly
[Wahl1988, Wahl2002](@cite). The retained charges are then invariant under ``Z \to Z_0 - Z``, and
the multiplicity ratio is exactly one half, which the parameterization pins. ``R_T(A_0/2) = 1``
follows under the charge-resolved inversion and under the ratio of means, whose two averages run
over one set of nuclides. The mean of ratios cannot reproduce it: the per-charge ratios come in
reciprocal pairs of equal weight, and the arithmetic mean of ``x`` and ``1/x`` exceeds one unless
``x = 1``, which leaves ``R_T(A_0/2)`` below one, by ``1 \times 10^{-5}`` to
``2 \times 10^{-4}`` on the four shipped domains. The identity at the
symmetric split therefore no longer singles out the ratio of means, as it did when that was the
only order compared with the mean of ratios.

A run checks both conditions — the invariance of the charge set, and ``R_T`` from ``r_\nu = 1/2``
at ``A_0/2`` — and records them, so that an input violating them is reported rather than silently
absorbed. A tabulated charge distribution may be given in place of Wahl's model; its polarization
at ``A_0/2`` is taken as tabulated unless `zero_polarization_at_symmetry` forces it to zero, as
the published extraction did.

## Description by joined segments

The ratio extracted point by point is scattered, and for some datasets sparse, so it is
``r_\nu`` that is described by joined segments and ``R_T`` that follows by the exact relation above. The model is
continuous and piecewise-linear, written in the truncated-power basis

```math
f(x) = \beta_0 + \beta_1 (x - x_0) + \sum_k \gamma_k (x - \psi_k)_+,
```

in which continuity at every breakpoint holds identically and the fit at fixed breakpoints is one
weighted linear least-squares solve, yielding the parameter covariance [Muggeo2003](@cite).
Breakpoints are restricted to abscissae present in the data and searched exhaustively, which is
tractable over a fragment mass range and returns the global optimum. The number of segments is
chosen by the Bayesian information criterion [Schwarz1978](@cite), which prices both the extra
slope and the extra breakpoint; the criterion for every order examined is retained, so the choice
can be audited.

The coefficient covariance ``(X^\mathsf{T} W X)^{-1}`` is scaled by
``\max(1, \chi^2/\mathrm{dof})`` where the data quote uncertainties: it is inflated where the
residuals exceed what the quoted uncertainties predict, and never shrunk below the quoted scale.
Where no point quotes an uncertainty the weights are uniform and carry no scale, so the fit is
ordinary least squares and the covariance is scaled by ``\chi^2/\mathrm{dof}`` alone.

Each segment must hold `min_points_per_segment` data points and extend over at least
`min_segment_span` mass units from its first pivot to its last. At four points per segment on
consecutive mass numbers the two guards coincide; the span guard acts where abscissae repeat, as
in a combined curve, or where the point count is set lower.

A dataset is offered as a segmented curve of its own only if it provides a complete fragment pair
at a fraction `min_dataset_coverage` of the mass numbers of the fragmentation range. Below that
floor it is still read, diagnosed and pooled into the systematic trend, but no curve is fitted to
it alone: with pairs at few mass numbers the breakpoint search cannot place the minimum where the
data do not reach, and the curve it returns asserts structure between the measurements that a
consuming code could not tell from a measured feature.

Two constraints are enforced during the search because they are exact, not preferences: the fit
is pinned to one half at the symmetric split — for a dataset whose complete pairs begin above
it the identity has no abscissa to act on, and that curve is fitted unpinned — and it may not leave the interval ``(0, 1)``, outside which the
temperature ratio relation is undefined. Since a piecewise-linear function attains its extrema at
its pivots, the bound is tested exactly rather than sampled.

## One curve per dataset, and one systematic trend

The datasets of a fissioning nucleus can differ well beyond their quoted uncertainties, and where
they do, the temperature ratio can only be determined separately for each of them. A run therefore
produces one segmented curve per dataset that reaches the coverage floor, fitted to that dataset
alone.

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
each ratio by the yield of a different fragmentation. The normalization of ``Y`` cancels. The two fragments of a split are counted in one event, so pre-neutron
``Y(A) = Y(A_0 - A)`` exactly; by default each mass of a distribution measured on both wings takes
the mean of the two before the average, so that a difference between the wings — a backing loss
on one side of a double-energy measurement — enters once and alike. The published averages took
the distributions as measured.

This is the quantity the literature tabulates, and it is not the mean over the fragment mass
range. That mean weights every mass number equally, so the far-asymmetric tail, where the yield is
smaller by orders of magnitude, counts as much as the peak. Both are reported, under names that
distinguish them.

The value is unchanged by how the uncertainty is formed. The uncertainty propagates the fit
covariance,

```math
\sigma^2 = w^\mathsf{T} C\, w
         + \sum_{A_H} \left(\frac{R_T - \langle R_T \rangle}{\sum Y}\right)^2 \sigma_Y^2,
\qquad w = \frac{Y}{\sum Y}, \qquad C = D\, J \Sigma J^\mathsf{T} D,
```

with ``J \Sigma J^\mathsf{T}`` the covariance of the fitted ``r_\nu`` at the tabulated mass numbers
and ``D`` the diagonal matrix of the slope ``\partial R_T/\partial r_\nu`` of the inversion at
each of them. The
tabulated points of a fitted curve are functions of a few coefficients and are not independent.
The independent-points form, which replaces ``w^\mathsf{T} C w`` by
``\sum (Y/\sum Y)^2 \sigma_{R_T}^2``, is the approximation of the published tables; it is written
beside the propagated uncertainty and understates it by a factor of two to five for the curves of
the four shipped systems. The yield term is common to both, which is why a multiplicity dataset
quoting no uncertainties still yields an uncertain average. Correlations between the yields at
different mass numbers are neglected, the sources not reporting them.

## Combining datasets

Several measurements of one fissioning system are not merged into one. Each is fitted on its own,
and a further curve is fitted to their combination.

That combination is not a concatenation. At each mass number the available values are combined
with inverse-variance weights carrying an additional between-dataset variance ``\tau^2``,
estimated from their dispersion [DerSimonian1986](@cite),

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
