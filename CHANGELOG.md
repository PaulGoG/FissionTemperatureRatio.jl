# Changelog

Notable changes to FissionTemperatureRatio.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html). Before 1.0 a breaking change raises
the minor version and a patch release carries none.

## [Unreleased]

## [0.2.3] - 2026-10-02

No trend and no dataset curve changes. The run identifiers of ²³⁵U and ²³⁹Pu change in their
`excl` token, an exclusion having been added to each. Corrections after a review of 0.2.2 are
included; those of them that would break an accepted configuration or a written file are left for
0.3.0.

### Added

- `pooled_datasets(result)`, the labels of the datasets a systematic trend was combined from:
  those that form at least one fragment pair in the range and that the configuration does not
  exclude.
- `deviation_autocorrelation(result)`: the lag-one autocorrelation, along the mass axis, of the
  pooled datasets' deviations from the combined curve. The covariance of a trend takes the
  combined points as independent, so its written uncertainty, and that of its total average, is
  low by about `√((1 + ρ)/(1 − ρ))`: `ρ` is 0.71, 0.61, 0.70 and 0.44 for ²⁵²Cf, ²³⁵U, ²³⁹Pu and
  ²³³U, factors of 2.4, 2.0, 2.4 and 1.6. The written uncertainty is not enlarged. `ρ` and the
  factor are written under `[result.trend_uncertainty]` in `metadata.toml`, and `ρ` as a last
  column, `deviation_autocorrelation`, of `segmented_curves_<run identifier>.csv`. The method
  page states the consequence where it explains the low chi-squared of a trend.

### Changed

- Batenkov 2004 for ²³⁵U (EXFOR 41502005) and for ²³⁹Pu (41502006) are excluded by accession in
  the shipped configurations, with the reason. Each gives 21 masses on a grid of 4 u on which no
  mass has its complement, so neither formed a fragment pair or entered a pool before, and they
  are not paired against interpolated complements; the paper, AIP Conf. Proc. 769, 1003 (2005),
  doi:10.1063/1.1945175, presents them as preliminary, with statistical uncertainties only. The
  exclusion puts that on record in the configuration, the run metadata and the dataset
  diagnostics.
- A `[yield] exclude` without `subdirectory` has nothing to act on and nothing to be checked
  against. It is still read for its form, with a warning, but no longer kept in the
  configuration or written to the metadata, where an accession that names nothing could appear.
- On a clone without the measured input, the testsets that need it are reported as skipped in
  the test summary; they were passed over without a trace.
- The README and the method page say that the paper's segments were drawn by hand and that the
  paper did not use this package; the README had the segment count "chosen by the Bayesian
  information criterion as published".

### Fixed

- The `pooled` column of `dataset_diagnostics.csv` was true for a dataset that forms no fragment
  pair, which contributes nothing to a trend: four datasets of the shipped ²³⁵U and ²³⁹Pu runs.
  It is now true only for the datasets `pooled_datasets` lists.
- Two files of one directory carrying one accession, as a retrieval under another spelling of
  the author leaves them, were read and pooled as two measurements, and an exclusion by that
  accession applied to one of them. Such a directory is now refused, when the configuration is
  loaded and when the directory is read.
- An exclusion that names no dataset was refused, in a run started from a configuration built by
  hand, only after every dataset had been fitted. It is refused before anything is fitted, the
  exclusions under `[yield]` included.
- Where a pool quotes no uncertainty, the mean was weighted by the pooling factors and the
  dispersion about it was not. The dispersion is now `Σ f (r − r̄)²/(Σf − Σf²/Σf)`, the estimate
  of `τ²` for values of no quoted variance; equal factors give the sample variance as before. No
  mass of the four shipped trends is in this case.

## [0.2.2] - 2026-10-02

Validated on the ExforFissionData v0.2.4 retrieval, whose tables for the four systems are those
of v0.2.3 byte for byte.

This release broke the rule it states. Keying an exclusion on the accession changed the `excl`
token of the run identifier for the same exclusion, and with it the names of the run directory,
of the manifest, of the curve table and of the total-average file; and
`Configuration.excluded_datasets`, `Configuration.excluded_mass_yields` and the exclusions in the
run metadata became keyed by accession where they had been keyed by label. Both are breaking
changes and belonged in 0.3.0. The deprecated label form of an exclusion can still name the wrong
dataset, where a directory holds only another subentry of the same author and year; it goes in
0.3.0.

### Changed

- Zeynalov 2019 for ²³⁵U (EXFOR 41738002), from the paper and analysis already excluded for
  ²⁵²Cf, is kept out of the ²³⁵U pool by the shipped configuration, with the reason; its own curve
  is still fitted and written. Weighted with the Al-Adili 2016 `Y(A)` over `A_H` = 130 to 150,
  its light- and heavy-fragment multiplicities are 1.02 and 1.36 against 1.36 to 1.66 and 1.02 to
  1.17 in Fraser 1966, Nishio 1998, Al-Adili 2020 and Vorobyev 2010, with the pair sum, 2.39, on
  the scale of `ν̄`: the neutrons are divided wrongly between the wings, and `r_ν` is 0.571
  against 0.413 to 0.436. The ²³⁵U trend moves from 6 segments, `χ²/dof` = 0.19 and
  `⟨R_T⟩` = 1.0743 to 5 segments, 0.44 and 1.1230.

### Deprecated

- `{ dataset = "<label>", reason = … }` in an `exclude` list, for a dataset that has an EXFOR
  accession: it still resolves, with a warning naming the accession to write. The form remains
  for a tabulation that carries no accession. Its removal for archive datasets will come with a
  minor version.

### Fixed

- A pool of values that quote no uncertainty carried no pooling fraction in the chi-squared of
  the trend: 0.2.1 took the fraction out of a standard error, `s/√n_eff`, that did not hold it,
  where 0.2.0 had it once. The mean of values counting for `f` of a measurement each has the
  variance `s²/Σf`; the standard error is now `s/√(Σf)` and the fraction enters once, as in the
  quoted case. No mass of the four shipped trends is in this case.
- Identical unquoted values at one mass gave the combined point a zero uncertainty, which
  `fit_weights` refuses, and stopped the run. The point now quotes no uncertainty and takes the
  median weight.
- The between-dataset variance was estimated as `(Q − (k − 1))/C` with the fixed-effect weights
  `f/σ²`, for which the expectation of `Q` without such a variance is `Σf − Σ(f²/σ²)/Σ(f/σ²)`,
  equal to `k − 1` only where every `f` is one. `τ²` was low by about `σ²(1 − f)/f`, and
  truncated to zero where it should not have been, wherever an interpolated dataset is pooled.
  It is now estimated against that expectation, the method of moments for general weights of
  DerSimonian and Kacker, doi:10.1016/j.cct.2006.04.004. The four trends keep their segments and
  breakpoints. The combined `r_ν` moves at 38 of 49 masses for ²⁵²Cf (by at most 7 × 10⁻⁵), at
  31 of 41 for ²³⁹Pu (`A_H` = 123 to 160, by at most 0.003, its standard error by up to 22 %) and
  at 12 of 43 for ²³³U (`A_H` = 118 to 130, at most 4 × 10⁻⁴). The trend `R_T(A_H)` moves by less
  than 0.001 for ²⁵²Cf and ²³³U; for ²³⁹Pu by more than 0.002 over `A_H` = 122 to 131, at most
  0.026 at `A_H` = 126 (1.741 to 1.767), 0.4 of its uncertainty there. `⟨R_T⟩`: ²⁵²Cf 1.0879 to
  1.0878, ²³⁹Pu 1.0267 to 1.0268, ²³³U 1.0724 unchanged.
- An exclusion matched the display label, which gains an accession suffix as soon as a second
  dataset of the same author and year is retrieved; the exclusion then matched neither and the
  dataset returned to the pool with a warning. An exclusion now names its dataset by EXFOR
  accession, `{ accession = "41739002", reason = … }`, under `[multiplicity]` and `[yield]`
  alike, and one that names no dataset of its directory is an error when the configuration is
  loaded, as is one that names no dataset read when a run starts. The run identifier hashes the
  accessions, so the `excl` token of the ²⁵²Cf run changes.

## [0.2.1] - 2026-10-02

This release carries breaking changes, marked below, and by the rule above should have been
0.3.0.

### Added

- `[retrieval] min_package_version`, default `"0.2.3"`: the lowest ExforFissionData version whose
  retrieval records a run accepts. A record below it, or one stating no `[run] package_version`,
  is refused when the configuration is loaded, and input directories written by different
  versions are reported in a warning. ExforFissionData 0.2.3 put the ²³⁹Pu multiplicities of
  22650004 on the scale of their subentry and refuses the ²⁵²Cf yields 41425015 and 41425016.
  `check_retrieval_versions`, `Configuration.min_retrieval_version`, `min_package_version` in the
  run metadata.
- The manifest records the EXFOR accession of every dataset curve, as `accession` in its
  `[[segmented_curve]]` entry; the systematic trend has none. The existing keys keep their names,
  and the label still selects a curve. `curve_accessions(result)` gives the accession of every
  dataset by label, and `segmented_curves_<run identifier>.csv` and `dataset_diagnostics.csv`
  carry it as a column.

### Changed

- Zeynalov 2019 (EXFOR 41739002) is kept out of the ²⁵²Cf pool by the shipped configuration, with
  the reason; its own curve is still fitted and written. Its `r_ν` lies 0.03 to 0.08 above the
  consensus of the other fourteen sets over `A_H` = 130 to 160, from a deficit of light-fragment
  neutrons (`⟨ν_L⟩` = 1.71 against 2.01 to 2.12, pair sum 3.47 against `ν̄` = 3.76) that does not
  cancel in the ratio; the README gives the attribution. The ²⁵²Cf trend moves from 4 segments,
  `χ²/dof` = 1.88 and `⟨R_T⟩` = 1.0788 to 5 segments, 1.55 and 1.0879.
- **Breaking:** `segments.min_dataset_coverage` gated two measures and is replaced by two keys,
  each validated on load and each a token of the run identifier: `segments.min_pair_coverage`
  (`cov`), the fraction of the heavy masses of the range a `ν(A)` dataset must pair to offer a
  curve of its own, and `yield.min_yield_coverage` (`Ycov`, on a run that names a yield
  distribution), the share of the primary distribution's heavy-fragment yield a `Y(A)` must cover
  to be averaged over. Both default to 0.3. The retired key is refused with a message naming its
  replacements. `SegmentSettings.min_pair_coverage`, `Configuration.min_yield_coverage`.
- **Breaking:** the files of a dataset are named by the stem of its input file,
  `<accession>_<Author>_<year>`, in place of its label: `R_T_vs_A_H_segmented_22660005_K.Nishio_1998.csv`
  for what was `R_T_vs_A_H_segmented_K._Nishio_1998.csv`, and likewise the point-by-point tables,
  the pivots and the per-dataset figures. Labels, the `systematic_trend` label and its file names
  are unchanged; a consuming code reads the file names from the manifest.

- FissionFragmentsDomain.jl v0.2.3 in every environment, `[compat]` 0.2.3. The symmetrization of
  a mass yield distribution is its `symmetrized_yield(yields::MassYield, compound_mass)`, of the
  same definition as the `symmetrized_mass_yield` it replaces: the mean of two measured
  complements with `√(σ_A² + σ_{A₀−A}²)/2`, an unquoted uncertainty contributing nothing, an
  unpaired mass standing for its complement. No result changes.

### Fixed

- The pooling fraction of an interpolated dataset entered the chi-squared of the systematic trend
  twice. The combined curve carries the fraction `f` in its standard error, `σ/√f` for a single
  value and `(Σ f/(σ² + τ²))^(-1/2)` otherwise, and the segmented fit multiplied each residual by
  the measured fraction again, so a residual of the trend entered with `f²/σ²` where a
  measurement of a dataset enters with `f/σ²`. The trend's `wrss`, its reduced chi-squared, the
  criterion that selects the segment count and the covariance scale were low wherever an
  interpolated dataset was pooled. `fit_segments` now forms one weight `f/σ²` per point, `σ`
  being its uncertainty as one measurement, for the solve, for chi-squared and for the
  information matrix alike, with `Σ f` behind the degrees of freedom, and the trend is fitted to
  the combined values at the uncertainty of one measurement. The trend of a pool holding one
  interpolated dataset reproduces that dataset's own fit exactly; a pool of integer-mass
  datasets is unchanged. Selected segment counts are as before for the four shipped systems; the
  reduced chi-squared of the trend rises from 1.77 to 1.88 (²⁵²Cf, same pool), 0.186 to 0.192
  (²³⁵U), 0.42 to 0.49 (²³⁹Pu) and 1.67 to 1.71 (²³³U), and the ²³³U trend moves its second
  breakpoint from 140 to 141 and its `⟨R_T⟩` from 1.0740 to 1.0724.
- With the fraction in the information matrix, the coefficient covariance of an interpolated
  dataset's own curve is divided by `f`: its uncertainties grow by `1/√f`. Coefficients and
  chi-squared of those curves are unchanged.

- The run metadata could state the version of an earlier release, `[source] package_version`
  being read when the module is compiled and the project file not being a dependency of the
  compiled cache: a release that changed only the version kept the previous one. The project file
  is now declared with `include_dependency`.

### Removed

- `symmetrized_mass_yield`; use `FissionFragmentsDomain.symmetrized_yield`.

## [0.2.0] - 2026-10-01

### Added

- The fragmentation domain is FissionFragmentsDomain.jl's (v0.2.1), shared with the prompt
  emission codes that read `R_T(A_H)`: nuclides and systems, the shipped AME2020 mass table and
  Gilbert-Cameron shell corrections, Wahl's charge model, the fragmentation domain, the level density models, and the relation
  between `R_T` and `E*_H/TXE` with its inverse and slope. Every environment takes the package
  from its repository at the tag through `[sources]`, and every activation script sets
  `JULIA_PKG_USE_CLI_GIT`.
- `ratio_averaging = "charge_resolved"`, the exact inverse of a partition that gives every
  fragmentation its own `a_L/a_H`, with the relation reduced over the charge distribution rather
  than through an effective ratio.
- `level_density.mean_kinetic_energy_file`, a pre-neutron `⟨TKE⟩(A)` dataset from the EXFOR
  retrieval, weighting every fragmentation of the charge-resolved inversion by its mean total
  excitation `Q + E*_CN − ⟨TKE⟩(A_H)`. Validated on load: listed in a retrieval record beside it,
  reaching the heavy-mass range. Where both `A` and `A₀ − A` are measured the pair takes their
  mean, and a light-wing value alone stands for its heavy complement; interior gaps are
  interpolated, and masses beyond the measured span take the nearest measured value. All are listed
  in `metadata.toml` with the retrieval record, its parser revision and its SHA-1, and with the
  yield-weighted mean `⟨TKE⟩` and its offset from the energy standard of the system,
  `recommended_mean_total_kinetic_energy`, over every yield distribution read, with the heavy
  masses that distribution covers: a partial one biases the mean towards its own masses. The 233-U
  configuration names Geltenbort 1985. `MeanKineticEnergy`, `read_mean_kinetic_energy`,
  `mean_kinetic_energy_offset`.
- Reaction-code qualifiers from the retrieval records. `DERIV` and `SPA` datasets are used and
  flagged: in the log, in the new `qualifiers` and `flagged` columns of `dataset_diagnostics.csv`,
  and in `metadata.toml`. `RetrievalRecord`, `retrieval_record`, `retrieval_qualifiers`,
  `FLAGGED_QUALIFIERS`.
- `fragmentation.charge_distribution_file` takes `"wahl"`, the default, `"mean"`, or a tabulated
  `A dZ sigma_Z` read with FissionFragmentsDomain's reader; `zero_polarization_at_symmetry` sets
  `ΔZ(A₀/2) = 0` for the last two, as the published extraction did, and is refused with `"wahl"`.
- `level_density.deformed_branch`, default `true`: Gilbert-Cameron with its eq. (21) for deformed
  nuclei beside eq. (20); `false` applies eq. (20) throughout, the published Table 2 setting.
- `yield.symmetrize`, default `true`: the total average is taken over `Y(A)` with the pre-neutron
  identity `Y(A) = Y(A₀ − A)` imposed, each mass taking the mean of the two wings where both are
  measured and a mass measured on one wing alone standing for its complement;
  `symmetrized_mass_yield`. The published settings keep `false`. Over the distributions measured
  on both wings it moves `⟨R_T⟩` by −2.5 % to +1.2 %: over Göök's ²⁵²Cf distribution, which the
  ²⁵²Cf configuration takes, by −0.1 % to −0.3 %, and over Straede's ²³⁵U one by −0.2 % to −1.0 %.
- Where datasets are pooled into the systematic trend, a dataset whose masses the retrieval
  interpolated onto the integers is weighted by its measured points over its written rows,
  `mass_values_non_integer / rows_written` from its retrieval record: neighbouring interpolated
  rows share their bracketing points. `consensus(curves; weights)`, `pooling_weight`, a
  `pooling_weight` column in `dataset_diagnostics.csv`, `pooling_weights` in the run metadata.
- The segment count is selected with the sample size and the degrees of freedom counted in
  measurements: `fit_segments(...; measured)` takes the measured fraction of every point, an
  interpolated dataset's raw points over rows, and a combined point the weighted average of its
  contributors'. `SegmentedFit` gains `points` and `measured_points`; `dof` is real.
- The systematic trend is refitted with one and two segments more than selected; its total
  averages at those orders are the columns `R_T_one_more_segment` and `R_T_two_more_segments` of
  `total_average_R_T_<run identifier>.csv` and `segment_count_sensitivity` in the run metadata.
- `[yield] mass_yield_file`, one distribution in place of `subdirectory`: the shipped
  configurations take Y(A) and ⟨TKE⟩(A) from one primary experiment per system, named by EXFOR
  accession in the run metadata (`mass_yield_accessions`, `mean_kinetic_energy.accession`):
  Göök 2014 for ²⁵²Cf, Al-Adili 2016 for ²³⁵U, Wagemans 1984 for ²³⁹Pu and Geltenbort 1985 for
  ²³³U.
- The coverage of a yield distribution is measured in yield, against the primary distribution of
  the system: the primary's heavy-fragment yield at the masses the distribution holds, over its
  yield across the fragmentation range. A distribution below `min_dataset_coverage` is read but
  not averaged over: its average would describe those masses, not the fission yield.
  `mass_yield_coverage(yields, reference, heavy_masses)`, `ExtractionResult.mass_yield_coverage`,
  and `mass_yield_coverage` and `mass_yields_not_averaged`, with the reason, in the run metadata.
- `[yield] mass_yield_file` is the primary distribution and is required in `[yield]`: averaged
  over alone, or beside `[yield] subdirectory` the reference of every distribution's coverage. A
  directory enters the run identifier as one token, `Y`, hashed together with its reference and
  its exclusions.
- `[yield] exclude`, distributions of the directory kept out of every average, each with its
  reason; they are still read and reported.
- The pair-sum scale the retrieval records for a multiplicity read by its complement test is
  reported: `pair_sum_scale`, the `pair_sum_deviation` and `scale_consistent` columns of
  `dataset_diagnostics.csv`, `scale_inconsistent_datasets` in the run metadata, and a warning for
  a dataset off the scale of `ν̄`. Such a dataset is used, not corrected: a uniform scale cancels
  in `r_ν`.
- Every total average states its yield fraction, the share of the distribution's yield over the
  fragmentation range that falls at the curve's mass numbers: `yield_fraction`,
  `TotalAverage.yield_fraction`, `TotalAverage(curve, yields, heavy_masses)`, and a
  `yield_fraction` column of `total_average_R_T_<run identifier>.csv`.
- `segmented_curves_<run identifier>.csv`, one row per manifest curve keyed by its label: segments,
  pin, span, pairs, coverage, reduced chi-squared, imputed weights, range mean. `metadata.toml`
  names it under `[outputs]`.
- The manifest records `[domain]`: level density model and Gilbert-Cameron branch, averaging and
  excitation weighting, charges per mass, charge model, mass table, FissionFragmentsDomain version.
  `manifest_domain(result)`.
- `ExtractedCurve` carries `kind`, `pairs` and `coverage`.
- `write_results` refuses a run whose identifier would make a file name longer than 255 bytes.

- A run is one directory. `scripts/run.jl` writes the tables, the manifest and the provenance
  record to `data/sims/<system>/<run identifier>/` and the figures to
  `plots/<system>/<run identifier>/`. An existing run of the same identifier is moved aside as
  `<run identifier>#1`, `#2`, …, so no run is overwritten. The provenance record `metadata.toml`
  holds `run_metadata(result)`, the identifier, the configuration and script paths, a timestamp,
  the git commit through DrWatson's `tag!`, and the platform: Julia version, hostname, kernel,
  machine, CPU model and threads, Julia and BLAS threads, memory, and `versioninfo()`. A copy of the
  configuration sits beside it as `configuration.toml`.
- `segments.min_segment_span`, the smallest extent of a segment in mass units from its first pivot
  to its last, 1 to 50, default 3. It acts where abscissae repeat or the point count is set below
  four; otherwise it coincides with the point guard.
- `segments.min_dataset_coverage`, the fraction of the mass numbers of the fragmentation range at
  which a dataset must provide a complete pair to be offered as a segmented curve, 0 to 1, default
  0.3. A dataset below it is read, diagnosed and pooled into the trend, but no curve is fitted to
  it alone: the breakpoint search cannot place the minimum where the data do not reach, and the
  curve it returns would assert structure between the measurements. All four shipped
  configurations set both keys explicitly.
- `covariance(fit, xs)`, the covariance `J Σ Jᵀ` of the fitted values at the abscissae `xs`.
- `SegmentedCurve(label, fit, R_a)` builds a curve and carries `R_T_covariance`, the covariance of
  the tabulated temperature ratio, whose diagonal is the square of the tabulated
  `R_T_uncertainty`. `total_average(curve::SegmentedCurve, yields)` and `range_mean(curve)`
  propagate it; `TotalAverage` holds `value`, `uncertainty` and `uncertainty_independent_points`,
  and `ExtractionResult.total_average_R_T` maps curve label and yield label to one.
- `R_T_uncertainty_independent_points`, a column of `total_average_R_T.csv` beside the
  covariance-propagated `R_T_uncertainty`: the independent-points form of the published tables,
  `total_average(curve.R_T, yields)`, kept for comparison with them.
- `SymmetryDiagnostics`, with fields `charge_set_invariant`, `R_a_at_symmetric_split`,
  `pinned_curves` and `warnings`, and `ExtractionResult.dataset_outcomes`, stating for every
  dataset whether it offers a segmented curve and, if not, why. `DatasetDiagnostics` gains
  `coverage`, and `diagnose` takes the fragmentation range: `diagnose(data, ratio, A₀, A_H_range)`.
- Manifest fields. Each `[[segmented_curve]]` entry carries `breakpoints`, `pivot_A_H`,
  `pivot_r_nu`, `reduced_chi_squared`, `first_A_H`, `last_A_H`, `pairs`, `coverage`,
  `range_mean_R_T` and `total_average_R_T`; `[run]` carries `identifier`, `package_version` and
  `total_average_R_T_columns`.
- A DOI-keyed bibliography for the documentation through DocumenterCitations, cited inline on the
  method page and listed on a references page.
- `windows_apply_to_datasets`, extending the physics windows that place a breakpoint at the heavy
  magic fragment from the systematic-trend curve to every dataset.
- `heavy_mass_min`, the smallest heavy-fragment mass number of the fragmentation range. It
  defaults to the symmetric split, which is what the range always started at before. Pinning is
  refused unless the range does start there, since the fit is pinned at its first abscissa and a
  range starting higher would fix the ratio to one half where nothing says it is.
- `system_notation`, the typeset form of a fissioning system — `²⁵²Cf(sf)`, `²³³U(nth,f)` — for a
  figure label or a caption, kept separate from the path token `system_label` so that one name
  does not mean both.
- `RUN_IDENTIFIER_ABBREVIATIONS`, the one table a run-identifier token is abbreviated through,
  keyed by the key's dotted path in the configuration file: `level_density.ratio_averaging => avg`,
  `segments.min_dataset_coverage => mincov`. `run_identifier` refuses a key with no entry rather
  than inventing a token, and the test suite asserts that every key it uses has one.
- `docs/src/naming.md`, the vocabulary this package applies: quantities, identifier rules,
  configuration rules, file layout and column headers.
- `write_results(result, directory)`, which writes the tables and the manifest of a result into a
  directory of the caller's choosing and refuses one that already holds files, so a caller using
  the package as a library chooses where results go rather than inheriting the active project.
  `run_pipeline` writes nothing.
- A test that consumes a run the way a downstream code would — reading the manifest and the files
  it names, using no internal function and taking columns by position — and checks the contract
  the README documents: the system is identifiable without parsing a label, every segmented curve
  names a file that exists and parses, the temperature ratio is tabulated densely enough that
  interpolation is exact, and the identity at the symmetric split survives the round trip.
- A section of the README describing that contract, including the things a consumer can get
  wrong: matching on header text instead of column position, reading the segment pivots instead of
  the tabulated temperature ratio, and leaving the averaging order at the published default when
  the consuming code re-expands with unaveraged level density parameters.

- `scripts/` has its own environment carrying CairoMakie, and `scripts/run.jl` activates it as its
  first statement: `julia scripts/run.jl config/<system>.toml`. A pipeline run no longer depends
  on a plotting stack installed in the default environment.
- `InsufficientDataError`, thrown by `fit_segments` when valid arguments meet data that cannot
  support a fit. The pipeline treats it as an outcome for that dataset; an `ArgumentError` from
  the fit now propagates instead of being reported as a missing curve.
- Run provenance records `versioninfo()` and takes the commit from the package source tree instead
  of the active project, and `scripts/run.jl` copies the resolved `Manifest.toml` of the scripts
  environment into the run directory.
- One multiplicity-ratio and one temperature-ratio figure per dataset, showing its points, its
  segmented fit with the uncertainty band, and the systematic trend as a guide. The overview
  figures now carry the data and the trend only; sixteen fits and bands in one axis were not
  readable.

### Changed

- **The quoted results rest on inputs retrieved with ExforFissionData.jl v0.2.3**, which
  interpolates masses the archive gives as non-integer values onto the integers rather than
  rounding them, reads the uncertainties of the archive subentry, admits multiplicities coded
  without the `FRG` tag by a complement test on their data — for ²⁵²Cf Budtz-Jørgensen 1988,
  Göök 2014, Zeynalov 2011, Britt 1964 and Piksaykin 1977, for ²³⁹Pu Tsuchiya 2000 and Batenkov
  2004, for ²³⁵U Batenkov 2004 — and refuses the ²⁵²Cf `Y(A)` of Vorobiev 2001, which is not the
  inclusive yield. Under the published settings Apalin's ²³³U row of Table 1 moves from +0.60 %
  to +1.89 %, Fraser's from +1.09 % to −0.41 %, and the ²⁵²Cf Budtz-Jørgensen row from −0.36 % to
  −0.53 %; the ²³⁹Pu rows over Nishio 1995 `Y(A)` come back to within 0.7 %, where Apalin and
  Zamyatnin had deviated by 2 %, and are now asserted, the yield distribution inferred by
  agreement.
- **The default inversion is `"charge_resolved"`.** Along the systematic trends `R_T` moves by up
  to 3.6 × 10⁻³ from the ratio of means, at the doubly magic heavy fragment, and every `⟨R_T⟩` by
  −0.02 % to −0.06 %. The excitation weight moves ²³³U `⟨R_T⟩` by a further −0.11 % to −0.15 %.
  `"ratio_of_means"` stays selectable.
- **The charge distribution is Wahl's model by default**, the 1988 per-reaction parameters for the
  four shipped systems, which have `ΔZ(A₀/2) = 0` of their own. Under the ratio of means, `⟨R_T⟩`
  moves by +0.02 % to +0.06 % against the digitised tables, and for ²³³U by +0.26 % to +0.43 %
  against the mean values it took before.
- **The manifest is FissionFragmentsDomain's run record**, written by its
  `write_temperature_ratio_manifest`: `[system]`, `[run]` (ordinate, abscissa, columns),
  `[domain]`, and per curve its label, kind and two files. The per-curve fields it carried before
  are in `segmented_curves_<run identifier>.csv`, `total_average_R_T_<run identifier>.csv` and
  `metadata.toml`.
- **The systematic-trend label is `systematic_trend`**, `SYSTEMATIC_TREND_LABEL`, the token a
  consuming code looks it up by; the figures still read "Systematic trend".
- `total_average_R_T.csv` is `total_average_R_T_<run identifier>.csv`.
- The local `SegmentedCurve` is `ExtractedCurve`: the segmented fit of `r_ν`, the `R_T` derived
  from it and that curve's diagnostics. `SegmentedCurve` is FissionFragmentsDomain's tabulated
  curve, the one a consumer reads.
- `temperature_ratio` and `total_average` are methods of FissionFragmentsDomain's functions:
  `temperature_ratio(averaging, model, domain, r_ν::RatioCurve)`,
  `total_average(::RatioCurve, ::MassYield)` and `total_average(::ExtractedCurve, ::MassYield)`.
  Uncertainties propagate with the exact slope of the inversion.
- `SymmetryDiagnostics.R_a_at_symmetric_split` is `R_T_at_symmetric_split`, the inversion of
  `r_ν = 1/2` at `A₀/2`; the effective `R_a` is no longer part of a result.
- `level_density.mass_excess_file` defaults to `"ame2020"` and `level_density.shell_correction_file`
  to `"gilbert_cameron_1965"`, the tables FissionFragmentsDomain ships; the latter matches the
  RIPL copy read before in every tabulated value.
- `build_level_density_model(settings, masses)`; `build_mass_table`, `build_charge_model`.
- Run-identifier tokens: `chg`, `dZ0`, `mass`, `TKE`, `Ysym` added, and `sc` and `def` for a
  Gilbert-Cameron run; `dZ`, `sZ`, `dZfile` removed; `maxseg`, `minpts`, `minspan`, `windat`,
  `mincov` shortened to `nseg`, `npts`, `span`, `wdat`, `cov`, so that the file names repeating the
  identifier stay within 255 bytes for every shipped configuration under either model.
- `Configuration.system` is a `FissioningSystem`.

- The run identifier is DrWatson's `savename` over every configuration key that changes the
  result, tokens sorted and joined by `_`. The system is not a token, since it names the directory
  the identifier sits in, and `significant_digits` is not, since it changes rendering and not the
  number. A list or a path — the required windows, the exclusion list, the charge distribution
  file, the yield directory — enters as the first eight hexadecimal digits of the SHA-1 of its
  canonical spelling, paths relative to the data directory so that the same inputs staged on
  another machine give the same identifier, and is written in full into `metadata.toml` under
  `[identifier.hashed]`.
- The run identifier names the run directory and the manifest, `manifest_<run identifier>.toml`,
  and no other file: a table is named by its quantity, its abscissa and its label, and the
  summaries are `total_average_R_T.csv` and `dataset_diagnostics.csv`. A consuming code stages the
  whole run directory and selects the manifest by its prefix.
- The configuration loader refuses unknown sections and keys, naming them.
- An unquoted uncertainty is `missing` in `Multiplicity.σν`, `MassYield.σY` and `RatioCurve.σ`;
  zero denotes an exact value and occurs only at a pinned abscissa of a fitted curve. Readers store
  an absent, non-numeric or non-positive uncertainty as `missing`, and a multiplicity file may be
  `A nu` with no third column. `fit_weights` throws on a quoted zero. Where one fragment of a pair
  quotes no uncertainty its term is dropped and the propagated value, a lower bound, is carried as
  quoted, which reproduces the published values for such datasets; where neither does, the ratio's
  uncertainty is `missing` and the point takes the median weight of the quoted ones. Point-by-point
  tables write an unquoted uncertainty as an empty field, and the per-dataset figure prints "no
  quoted uncertainties" in place of a reduced chi-squared for a dataset quoting none.
- The uncertainty of the total average of a segmented curve propagates the fit covariance: the
  ratio term is `wᵀ C w`, with `w` the normalized yield weights, in place of a sum of squares over
  tabulated points that are not independent. The value is unchanged. **This moves the
  uncertainties**, by a factor of two to five for the shipped systems: ²³³U Nishio 1998 over Surin
  1972 is `1.1822 ± 0.0050`, against `± 0.0020` in the independent-points form and a published
  `1.1861 ± 0.0021`; Apalin 1965 is `1.0541 ± 0.0089`, against `± 0.0022` and a published
  `1.0346 ± 0.0043`. The range mean is propagated the same way.
- Figures use the standard layout — a 900 × 600 canvas per panel, 1200 wide above eight legend
  entries, 26 pt type, 3 pt data lines, 14 pt markers with a darker edge — in place of the
  journal-column sizing. Error bars are drawn without caps. `plot_ratio` gains `fit_label`,
  `yticks` and `annotation_corner`.

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
  `dataset_diagnostics.csv`. The manifest lists `[[segmented_curve]]` entries and declares the
  ordinate, the abscissa and the column names.
- Caught exceptions are bound to `exception` rather than `err`, which the naming convention
  forbids as an abbreviation and which the two companion packages do not use either.
- `Printf` is a dependency, for writing a value and its uncertainty to the same precision.

### Fixed

- The covariance was scaled by χ²/dof even below one, shrinking the band below what the quoted
  uncertainties support. The coefficient covariance `(XᵀWX)⁻¹` is now scaled by `max(1, χ²/dof)`
  where the data quote uncertainties, and by `χ²/dof` alone where no point quotes one, since
  uniform weights carry no scale.
- Two datasets of one author and year — two subentries of one measurement — shared a label, and a
  label selects a curve downstream. Each now carries its archive identifier in parentheses. An
  identifier of nine digits, a pointer within a subentry, is stripped from the label like one of
  eight.
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
- Two yield distributions of one author and year shared a label, and the total averages are keyed
  by label, so the second replaced the first: over the 252-Cf distributions of Barreau 1985 the
  far-wing subentry 23717005 stood in for the full 23717003. Each now carries its archive
  identifier in parentheses, as the multiplicity datasets do.
- The published total averages are now asserted by the test suite, where the input data is
  present, instead of being compared by the documentation script only.
- The charge distribution tables are cited to *At. Data Nucl. Data Tables* **39**, 1 (1988); the
  volume was given as 38. The liquid-drop coefficients are attributed to Pearson, *Hyperfine
  Interact.* **132**, 59 (2001), which is where the 2005 level density paper takes them from, and
  literature references carry DOIs.
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

### Removed

- The local mass table reader, level density models, charge distribution, fragmentation domain,
  `average_over_charge`, `level_density_ratio`, `RatioAveraging` and its subtypes, `MassYield`,
  `read_mass_yield`, `mass_yield`, `REACTIONS`, `CHANNEL_REACTION`, `element_symbol`,
  `system_label`, `system_notation` and `TREND_LABEL`: FissionFragmentsDomain's are used instead,
  under the same names where they exist.
- The forced `ΔZ(A₀/2) = 0`; the charge model provides it, and `zero_polarization_at_symmetry`
  imposes it on a table.
- `fragmentation.fallback_charge_polarization` and `fallback_charge_dispersion`.

- `weighted_mean`. The summaries of a segmented curve, `range_mean` and `total_average`,
  propagate its covariance instead; an inverse-variance mean over correlated points had no
  meaning.
- `CITATION.cff`. The README names the author and cites the method paper by its DOI, and the
  repository by its URL.

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

[unreleased]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/compare/v0.2.3...HEAD
[0.2.3]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/PaulGoG/FissionTemperatureRatio.jl/releases/tag/v0.1.0
