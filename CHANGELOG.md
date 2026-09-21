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
- `windows_apply_to_datasets`, extending the physics windows that place a breakpoint at the heavy
  magic fragment from the systematic-trend curve to every dataset.
- `heavy_mass_min`, the smallest heavy-fragment mass number of the fragmentation range. It
  defaults to the symmetric split, which is what the range always started at before. Pinning is
  refused unless the range does start there, since the fit is pinned at its first abscissa and a
  range starting higher would fix the ratio to one half where nothing says it is.
- `system_notation`, the typeset form of a fissioning system — `²⁵²Cf(sf)`, `²³³U(nth,f)` — for a
  figure label or a caption, kept separate from the path token `system_label` so that one name
  does not mean both.
- `RUN_IDENTIFIER_ABBREVIATIONS`, the one table a run-identifier token is abbreviated through.
  `run_identifier` refuses a key with no entry rather than inventing a token, and the test suite
  asserts that every key it uses has one. `parsimony` now appears in the identifier, which it did
  not, although it changes the result.
- `docs/src/naming.md`, the vocabulary this package applies: quantities, identifier rules,
  configuration rules, file layout and column headers.
- The output root of a run is configurable, so a caller using the package as a library chooses
  where results go rather than inheriting the active project.
- A test that consumes a run the way a downstream code would — reading the manifest and the files
  it names, using no internal function and taking columns by position — and checks the contract
  the README documents: the system is identifiable without parsing a label, every segmented curve
  names a file that exists and parses, the temperature ratio is tabulated densely enough that
  interpolation is exact, and the identity at the symmetric split survives the round trip.
- A section of the README describing that contract, including the things a consumer can get
  wrong: matching on header text instead of column position, reading the segment pivots instead of
  the tabulated temperature ratio, and leaving the averaging order at the published default when
  the consuming code re-expands with unaveraged level density parameters.

### Changed

One name per quantity, from the configuration key to the column header. The package now shares its
vocabulary with the retrieval that supplies its input and the emission model that consumes its
output, so a name learned in one is the name in the others. This is a breaking change to the
configuration schema, to the exported names and to the layout of `data/`; no deprecated aliases
are provided, because the schema changed incompatibly regardless and a package that answers to old
names while refusing old keys is worse to debug than a clean break.

- A fissioning system is declared by its entrance channel — `sf`, `nth`, `nres`, `nfast` — rather
  than by a reaction code. The reaction code follows from the channel and is no longer a key: it
  could only be redundant or wrong, and `n,f` alone cannot separate a thermal run of a system from
  a resonance run of the same system, which the system identifier must. System identifiers are
  therefore `Cf252_sf`, `U233_nth`, `U235_nth`, `Pu239_nth`, and configuration files are named for
  the system they run.
- `output.digits` is now `output.significant_digits` and means significant figures rather than
  decimal places. **This moves the numbers**, by at most one part in 10⁵ on a ratio and by rather
  more on an uncertainty, which gains precision: an uncertainty of `5.49433e-4` was written
  `0.000549`, to three significant figures, and is now written in full.
- Configuration keys: `A_H_max` → `heavy_mass_max`; `fallback_rms` → `fallback_charge_dispersion`;
  `level_density.prescription` → `level_density.model`; `multiplicity.directory` and
  `yield.directory` → `subdirectory`, since `subdirectory` names a folder under a root and
  `directory` only ever names a root; `exclude`'s `set` → `dataset`.
- Types: `ChargeDistributionData` → `ChargeDistribution`, `MultiplicityData` → `Multiplicity`,
  `YieldData` → `MassYield`, `DataSetDiagnostics` → `DatasetDiagnostics`, `Parameterization` →
  `SegmentedCurve`, `PipelineResult` → `ExtractionResult`, `LevelDensityPrescription` →
  `LevelDensityModel`.
- Readers are named for what they return: `read_mass_excess` → `read_mass_excess_table`,
  `read_shell_corrections` → `read_shell_correction_table`, `read_yield` → `read_mass_yield`.
  `build_prescription` → `build_level_density_model`, `case_label` → `system_label`.
- `rms` is retired throughout. The quantity is the Gaussian dispersion of the isobaric charge
  distribution: `σ_Z` in code, `sigma_Z` in files and headers.
- `data/` is laid out as `reference/` for the system-independent evaluations and one directory per
  system, with one subdirectory per measured quantity — the layout the retrieval writes. File
  extensions no longer carry a nuclide.
- Column headers name their quantity and their uncertainty: `A nu nu_uncertainty`,
  `A_H,R_T,R_T_uncertainty`. No column is called `value`. Readers continue to take columns by
  position, which is what makes every header rename a no-op for code, and each reader's docstring
  now says so.
- Result files state the quantity and the abscissa: `R_T_parameterized_…` →
  `R_T_vs_A_H_segmented_…`, `segments_…` → `r_nu_vs_A_H_pivots_…`, `diagnostics_…` →
  `dataset_diagnostics_…`. The manifest lists `[[segmented_curve]]` entries and declares the
  ordinate, the abscissa and the column names.
- Caught exceptions are bound to `exception` rather than `err`, which the naming convention
  forbids as an abbreviation and which the two companion packages do not use either.
- `Printf` is a dependency, for writing a value and its uncertainty to the same precision.

### Fixed

- **The pin `r_ν = 1/2` was applied at the first abscissa of each dataset, not at the symmetric
  split.** `r_ν = 1/2` is an identity at `A₀/2` only, but a dataset whose first complete fragment
  pair lies above it was pinned to one half there: 235-U Nishio 1998 at `A_H = 126`, where the
  measured ratio is 0.28, 252-Cf Britt 1964 at 136, and eleven others across the four shipped
  systems. **This moves the numbers** of those thirteen curves — `R_T(A_H)` near the first
  abscissa by tens of per cent, `⟨R_T⟩` by up to 31 % for the sparsest sets — and leaves every
  curve that starts at the symmetric split, and every systematic-trend curve, unchanged. Such
  datasets are now fitted unpinned, the run reports which, and each manifest entry carries
  `pinned_at_symmetric_split`. Against Table 1 of the paper the 239-Pu Nishio row moves from
  2.0 % to 0.09 %, 235-U Nishio from 0.78 % and 1.05 % to 0.56 % and 1.00 %, and 233-U Fraser from
  0.27 % to 1.09 %.
- The published total averages are now asserted by the test suite, where the input data is
  present, instead of being compared by the documentation script only.
- Tick labels of adjacent panels ran together in the method-chain figure. The row gaps were being
  set by index before the legend was added, which renumbers them, so the gap that was tightened
  was not the one intended.
- The in-axis annotation of the temperature-ratio figure was built as a plain string, so it
  printed `⟨R_T⟩` with a literal underscore instead of a subscript. It is typeset now, one
  `LaTeXString` per line, since a line is one expression and the engine has no line break.
  `plot_ratio` accordingly accepts a list of lines for `annotation` as well as a single string.
- That annotation sat at the lower right, which for a temperature ratio is where the heavy wing
  descends: the curves ran straight through the text. It is anchored at the upper right, the one
  corner the quantity leaves empty in every system.
- The annotation printed a value and its uncertainty at different precisions — `1.18 ± 0.006` —
  because rounding drops a trailing zero. Both are written to the same three decimal places.

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
