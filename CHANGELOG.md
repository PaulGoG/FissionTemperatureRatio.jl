# Changelog

Notable changes to FissionTemperatureRatio.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `parsimony`, multiplying the penalty the selection criterion charges per parameter, as the knob
  for biasing the number of segments downwards. One is the criterion as published. A flat offset
  was tried first and rejected: it scales neither with sample size nor with how many parameters a
  segment costs, and on this data the criterion prefers five segments to four by a margin already
  counted as very strong evidence, so no conventional offset moves it.
- `windows_apply_to_data_sets`, extending the physics windows that place a breakpoint at the heavy
  magic fragment from the systematic-trend curve to every data set.

- The output root of a run is configurable, so a caller using the package as a library chooses
  where results go rather than inheriting the active project.
- A test that consumes a run the way a downstream code would — reading the manifest and the files
  it names, using no internal function — and checks the contract the README documents: the system
  is identifiable without parsing a label, every parameterization names a file that exists and
  parses, the temperature ratio is tabulated densely enough that interpolation is exact, and the
  identity at the symmetric split survives the round trip.
- A section of the README describing that contract, including the two things a consumer can get
  wrong: reading the segment pivots instead of the tabulated temperature ratio, and leaving the
  averaging order at the published default when the consuming code re-expands with unaveraged
  level density parameters.

## [0.1.0] - 2026-09-13

First release. The package extracts the temperature ratio of complementary fully accelerated
fission fragments from experimental prompt neutron multiplicity data, parameterizes it, and
reproduces the total averages published for 233-U(n,f), 235-U(n,f) and 252-Cf(sf) to better than
one per cent.

### Added

- A fissioning system is declared as target, reaction and incident energy; the fissioning nucleus
  and the case label are derived from them. A label can no longer contradict the nuclide it names,
  and two incident energies of one target are two systems, distinguished in the run identifier.
- Combination of several measurements of one system by inverse-variance weighting with a
  between-set variance term, mass number by mass number, in place of concatenating them. The sets
  disagree by ten to twenty times their quoted uncertainties, so a concatenation weighted by those
  uncertainties is decided by whichever author quoted the smallest ones. The reduced chi-squared of
  the trend curve falls from about 16 to about 2 for 252-Cf, and now describes the fit rather than
  the disagreement. The combined curve is written out.
- Structural diagnostics for every data set — usable pairs and their span, points outside the
  physical range, departure from one half at the symmetric split, complementary multiplicities
  against the total — written whether or not the set was used. Sets are kept out of the
  combination only when the configuration names them and gives a reason; an excluded set is still
  read, fitted, written and diagnosed.
- A run manifest naming the system, every parameterization, and the file to read for each, for
  codes that consume these curves rather than read them.
- A propagated uncertainty column in the segment files.

- Fragment mass yield input and the total average `⟨R_T⟩ = Σ Y(A_H) R_T(A_H) / Σ Y(A_H)`, reported
  for every combination of parameterization and yield distribution. This is the quantity the
  literature tabulates; the mean over the fragment mass range, which the run also reports, weights
  every mass number equally and is dominated by the far-asymmetric tail. The normalization of `Y`
  cancels. Both the ratio and the yield propagate into the uncertainty, which is what gives a
  multiplicity set quoting no uncertainties a finite one.
- Reproduction of the published total averages. Against Tables 1 and 2 of Eur. Phys. J. A **60**,
  190 (2024), from independently retrieved archive data: 233-U(n,f) to 0.60 %, 252-Cf(sf) to
  0.55 %, 235-U(n,f) to 1.05 %, the Gilbert-Cameron variant to 0.21 %. 239-Pu(n,f) deviates by up
  to 2.1 %, the published table having averaged over a yield distribution it does not name and a
  second that is calculated rather than measured.
- Experimental input sourced with `ExforFissionData.jl` in place of hand-assembled files. Every
  data file carries its EXFOR DatasetID and each directory holds the retrieval run record; the
  readers ignore anything that is not a `.dat`, so the record sits beside the data.
- The plotting stack as a weak dependency. `publication_theme`, `plot_multiplicities`,
  `plot_ratio`, `save_figure` and `write_figures` are declared by the package and implemented in
  an extension that loads with CairoMakie, so a package whose output is tabulated data no longer
  costs a graphics stack on `using`: about a second, with no Makie loaded. A run without it writes
  every table and its metadata, and says that figures were skipped.
- `CITATION.cff`, a version-bounded `formatter/` environment, a `check.jl` pre-commit gate, a
  JuliaFormatter CI workflow, and documentation deployment.
- Extraction of the temperature ratio `R_T = T_L/T_H` of complementary fully accelerated
  fragments from experimental prompt neutron multiplicity `ν(A)`, with no prompt emission
  calculation entering, following Eur. Phys. J. A **60**, 190 (2024).
- Parameterization of the multiplicity ratio `r_ν = ν_H/(ν_L + ν_H)` by joined straight segments.
  Breakpoints are restricted to abscissae present in the data and searched exhaustively, so the
  criterion attains its global optimum; the number of segments is chosen by the Bayesian
  information criterion, which prices the breakpoints as the fitted parameters they are.
- One parameterization per `ν(A)` data set, and separately a systematic-trend curve fitted
  through the pooled data with the minimum forced into `A_H` 128-132. The paper determines
  `R_T(A_H)` per data set where the sets disagree, and supplies the trend curve for the cases
  where the data cannot resolve the shape.
- Two level density prescriptions behind one interface, each carrying its own tabulated data:
  back-shifted Fermi gas and Gilbert-Cameron, selected by configuration.
- Both orders of reducing the level density parameter ratio over the isobaric charge
  distribution, `⟨a_L⟩/⟨a_H⟩` and `⟨a_L/a_H⟩`. The first is the default: it alone returns exactly
  one at the symmetric split, where the two fragments are the same nuclide.
- Validated TOML configurations for 233-U(n,f), 235-U(n,f), 239-Pu(n,f) and 252-Cf(sf). The
  parser enforces the types, enumerated choices and numerical bounds its comments document, and
  fails naming the offending key.
- Run provenance: a run identifier, the git commit, the resolved versions of the direct
  dependencies, and a hardware fingerprint, written with every result. Results are written
  through `safesave`, so no previous run is overwritten.
- Publication-width CairoMakie figures at journal column size.
- Uncertainty propagated through the segment model from the parameter covariance, scaled by the
  reduced chi-squared, so the band widens away from the data and vanishes at a pinned abscissa.
- Test suite with Aqua, JET restricted to the package's own frames, and ExplicitImports;
  benchmarks for the breakpoint search and the level density parameter ratio.

### Fixed

Corrections relative to the prototype this package replaces, preserved on the `legacy/prototype`
branch:

- The level density parameter ratio was formed per charge and then averaged, `⟨a_L/a_H⟩`, where
  the published method uses `⟨a_L⟩/⟨a_H⟩`. Both are now available and the published order is the
  default; the two differ by Jensen's inequality and coincide only for a degenerate charge
  distribution.
- The segment-selection objective was the mean residual sum of squares over segments, with no
  penalty on breakpoint placement, minimised from a magic initial score. It preferred more
  segments by construction.
- The uncertainty band carried the slope standard error alone, and was zeroed before the second
  fit, so it was not the propagated uncertainty of the segment model. Reported average
  uncertainties used `sqrt(Σσ²)/N`, which is not the uncertainty of a weighted mean.
- `ν(A)` points were deleted by undocumented thresholds on the ratio of the third column to the
  second, which guessed whether that column was an uncertainty. Experimental data is no longer
  discarded on the basis of its value; a point quoting no uncertainty is given the median weight
  and recorded as imputed.
- The multiplicity ratio was clamped to exactly one where it exceeded one, which is the singular
  point of the transformation to `R_T`. Such points are now excluded, and a candidate fit leaving
  the physical interval is rejected outright rather than repaired.
- The fragmentation domain sorted the mass and charge columns independently and then looked
  values up by pair, which can mispair rows and throws on duplicates. Rows are sorted as rows.
- The guard suppressing the duplicate push at the symmetric split compared the light fragment
  mass against the first mass of the range rather than against the heavy fragment mass. The two
  agree only when the range starts at `A₀/2`; the guard is now stated as `A_L != A_H`.
- Parameters were hardwired as constants in the source: target nucleus, `A₀`, `Z₀`, charges per
  mass number, segment count and bounds, rounding. They are configuration.
- Paths were relative to a `cd` at load, output directories were created on a single existence
  check so a partially present tree failed, and results were overwritten in place with no run
  identifier, commit or hardware record.

[unreleased]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/releases/tag/v0.1.0
