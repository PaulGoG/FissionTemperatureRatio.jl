```@meta
CurrentModule = FissionTemperatureRatio
```

# Naming

A quantity has exactly one name, and that name is used everywhere it appears: Julia identifier,
struct field, configuration key, directory name, file name, column header, figure label, prose.
The spelling changes register between contexts — `heavy_mass_max` in a configuration, `A_H` in a
path — but the *quantity* does not.

A name says what the thing is: `R_T_vs_A_H`, not `R_T_A_H`; no `Data` suffix; no column called
`value`.

This page gives the vocabulary as this package applies it. The retrieval that supplies its
experimental input and the prompt emission codes that consume its output share the quantity table
and the file layout, so a name learned here is the name there.

## Two registers

| Register | Spelling | Where |
| :--- | :--- | :--- |
| word | `heavy_mass_max`, `multiplicity` | configuration keys and values, prose |
| symbol | `A_H`, `nu` | directory names, file names, column headers, figure labels, Julia identifiers of the equations |

A configuration is edited by hand and has to explain itself, so it spells a quantity out. A path
and a header are read at a glance and are the field's own nomenclature, so they carry the symbol.
This is why [`FragmentationSettings`](@ref) has a field `heavy_mass_max`, mirroring the key it was
read from, while a ratio curve has a field `A_H`, mirroring the equations.

## The quantities

| Symbol | ASCII token | Julia | Quantity |
| :--- | :--- | :--- | :--- |
| *A* | `A` | `A` | pre-neutron fragment mass |
| *A*_H, *A*_L | `A_H`, `A_L` | `A_H`, `A_L` | heavy, light fragment mass |
| *A*₀, *Z*₀ | `A_0`, `Z_0` | `A₀`, `Z₀` | fissioning nucleus mass, charge |
| *Z* | `Z` | `Z` | fragment charge |
| *Z*_p | `Z_p` | `Zₚ` | most probable charge |
| Δ*Z* | `dZ` | `ΔZ` | charge polarization |
| σ_*Z* | `sigma_Z` | `σ_Z` | charge dispersion |
| *p*(*Z*,*A*) | `p_Z` | `p` | isobaric charge probability |
| Δ | `mass_excess` | `Δ` | mass excess |
| *S*(*N*), *S*(*Z*) | `S_N`, `S_Z` | `S_N`, `S_Z` | shell correction |
| δ*W* | `dW` | `δW` | shell correction plus pairing, as the Fermi-gas systematic uses it |
| *a* | `a` | `a` | level density parameter |
| ρ_*Z* | `rho` | `ρ` | level density parameter ratio of one fragmentation, a_L/a_H |
| *Q* | `Q` | `Q` | energy released in one fragmentation |
| TKE | `TKE` | `TKE` | pre-neutron total kinetic energy |
| TXE | `TXE` | `TXE` | total excitation energy, Q + E*_CN − TKE |
| ν(*A*) | `nu` | `ν` | prompt neutron multiplicity per fragment |
| *Y*(*A*) | `Y` | `Y` | pre-neutron fragment mass yield |
| *r*_ν | `r_nu` | `r_ν` | multiplicity ratio, ν_H/(ν_L + ν_H) |
| *R*_a | `R_a` | `R_a` | effective level density parameter ratio, ⟨a_L⟩/⟨a_H⟩ or ⟨a_L/a_H⟩ |
| *R*_T | `R_T` | `R_T` | temperature ratio, T_L/T_H |

`rms` is not used: the quantity is the Gaussian dispersion of the isobaric charge distribution,
σ_Z.

### Uncertainty

One spelling, one position: **`<quantity>_uncertainty`**, immediately after the quantity it
belongs to. `A nu nu_uncertainty`, never `A nu errnu`, and never a bare `uncertainty` column that
leaves the reader to work out which quantity it belongs to. In Julia the symbol is used instead —
`σν`, `σY` — because there the quantity it belongs to is not in doubt.

## Julia identifiers

| Kind | Case | Example |
| :--- | :--- | :--- |
| module, type, abstract type | `CamelCase` | [`ExtractedCurve`](@ref), [`SegmentedFit`](@ref) |
| function, macro | `lowercase_snake_case` | [`run_pipeline`](@ref), [`read_mean_kinetic_energy`](@ref) |
| variable, field, keyword | `lowercase_snake_case` | `charges_per_mass`, `significant_digits` |
| constant | `SCREAMING_SNAKE_CASE` | [`FLAGGED_QUALIFIERS`](@ref) |

Unicode is used where an identifier *is* a symbol of the governing equations — `ΔZ`, `σ_Z`, `ν`,
`Zₚ`, `A₀`, `R_T` — and only there, so that the implementation reads like the formulation. The
exported API stays ASCII-typeable: every public function can be called, and every keyword argument
of one given, without a compose key. That is why [`plot_multiplicities`](@ref) takes `A_0` while
the field it is passed is `A₀`.

Verb prefixes are a small closed set — `read_`, `write_`, `build_`, `load_`, `run_`, `is_` — and a
function that fits none of them is named for what it returns. `read_<thing>` must return the type
named `<Thing>`, which is the rule that decides `read_mass_excess_table` against
`read_mass_excess`: the type is `MassExcessTable`, so the reader carries `table` too.
`load_<thing>` reads *and* validates, so the configuration entry point is
[`load_configuration`](@ref) and not `read_configuration`.

Type names carry no `Data` suffix: a type holding multiplicities is [`Multiplicity`](@ref), one
holding mass yields is `MassYield`, one holding the charge polarization and dispersion is
`ChargeDistribution`. `Table` survives only where the thing genuinely is a lookup table keyed by
nucleon number. No adjective stands in for a noun, and no abstract noun stands in for a
description: the piecewise-linear fit of one measurement is a [`SegmentedFit`](@ref), and the
temperature ratio extracted through it an [`ExtractedCurve`](@ref). The tabulated curve a
consuming code reads is FissionFragmentsDomain's `SegmentedCurve`; the two names are distinct
because the two things are.

The fragmentation domain, the nuclides and systems, the mass table, the charge models, the level
density models, the relation between `R_T` and `E*_H/TXE`, the mass yield and the run record are
FissionFragmentsDomain's, and keep its names here; this package adds methods to its functions,
never a second definition under the same name.

The result type of a package is `<Verb>Result`, **one per package**. Here it is
[`ExtractionResult`](@ref), named for what this package does: it extracts a temperature ratio.
`RunResult` would collide, on `using` both, with the same name in a code that consumes these
curves.

Abbreviations are permitted only where the field itself uses them — `TKE`, `TXE`, `BSFG`, `GC`,
`EXFOR`, `AME`, `RIPL`, `sf`, `BIC` — and never as `err`, `param`, `calc`, `val`, `idx`, `num`,
`cfg`, `pts`. The one place short forms are unavoidable is the readable tokens of a run
identifier, seven of them; see below.

## Configuration files

A configuration is named `<system>.toml`, the token that names the data directories it reads and
the output directory it writes. A system carrying more than one configuration would name them
`<system>_<case>.toml`; none does. The test suite asserts the name against the system each file
declares, so a file cannot drift from what it runs.

Sections are `[system]`, `[fragmentation]`, `[level_density]`, `[multiplicity]`, `[yield]`,
`[segments]`, `[retrieval]` and `[output]`, and **a key is never prefixed with its section's
name**: `[level_density] model`, not `level_density_model`. Every key, with its type, bounds and
default, is listed under [Configuration](configuration.md).

A key's comment states three things and nothing else: what the key is, in one short line; the
enumerated choices where it is constrained; and the bounds and the unit. No physics narrative, no
run rationale — those belong in the documentation. The loader enforces exactly what the comment
declares, so a run cannot start from a configuration it cannot honour.

A key that can only be redundant or wrong is not in the file. The reaction code follows from the
entrance channel, so it is not written down; the fissioning nucleus follows from the target and
the channel, so it is not written down either.

`subdirectory` names a folder under a root; `directory` only ever names a root path. So a
configuration says `[multiplicity] subdirectory = "U233_nth/nu_vs_A"`, under the data root, and
the root itself is an argument to [`load_configuration`](@ref).

`significant_digits` is significant figures. It is not decimal places, and the distinction is not
cosmetic: a ratio near unity and an uncertainty near `10⁻⁵` appear in the same table, and a fixed
number of decimals cannot carry both to the same meaning.

## The system identifier

One token for a fissioning system, used as a directory name, a configuration file name and an
output directory:

```
<ElementSymbol><A>_<entrance channel>
```

with the channel spelled as the field spells it and not as a reaction code writes it — `Cf252_sf`,
`U235_nth`, `U235_nres`, `U233_nth`, `Pu239_nth`. `0f` is a reaction code and belongs in the run
record. The channel is also what distinguishes a thermal from a resonance run of one system, which
`n,f` alone cannot.

FissionFragmentsDomain's `system_label` returns that token. The typeset form for a figure or a
caption — `²⁵²Cf(sf)`, `²³³U(nth,f)` — is its `system_notation`, a separate function: one name may
not mean both.

## Directories, files and headers

```
data/
└── <system>/                         Cf252_sf, U233_nth, U235_nth, Pu239_nth
    ├── charge_distribution_vs_A.dat  optional: a tabulated charge distribution
    ├── nu_vs_A/
    │   ├── retrieval.toml            the run record of the retrieval
    │   └── <accession>_<Author>_<year>.dat
    ├── Y_vs_A/
    │   └── <accession>_<Author>_<year>.dat
    └── TKE_vs_A/
        └── <accession>_<Author>_<year>.dat
```

The atomic mass evaluation and the Gilbert-Cameron shell corrections are the ones
FissionFragmentsDomain ships, named `"ame2020"` and `"gilbert_cameron_1965"` in a configuration.

One directory per system, one subdirectory per measured quantity, one file per measurement. The
directory carries the quantity and the file carries the provenance, which is what makes the
archive accession number recoverable from a figure legend.

A tabular file is named `<quantity>_vs_<abscissa>[_<abscissa>]`. Extensions state a format and
nothing else: `.dat` whitespace-separated, `.csv` comma-separated, `.toml` metadata. A file
extension never carries a nuclide; the nuclide is in the path.

One header line, ASCII, abscissae first, then the ordinate, then its uncertainty:

```
A nu nu_uncertainty
A Y Y_uncertainty
A TKE TKE_uncertainty
A dZ sigma_Z
A_H,R_T,R_T_uncertainty
```

No column is named `value`: a column named for its quantity is readable without the file name,
which matters most for the files that cross a package boundary.

Readers take columns **by position**, not by header text. A header rename then changes nothing,
upstream or downstream. Every reader in this package says so in its docstring, and none may be
changed to a lookup by name.

## Output files and run identifiers

A run is one directory, and the directory carries the run identifier. Tables go under the system
in `data/sims/`, figures under the system in `plots/`:

```
data/sims/<system>/<run>/
├── manifest_<run>.toml                        what the run produced, for a consuming code
├── r_nu_vs_A_H_<dataset>.csv
├── R_T_vs_A_H_<dataset>.csv
├── r_nu_vs_A_H_consensus_systematic_trend.csv
├── r_nu_vs_A_H_segmented_<dataset>.csv        …_systematic_trend.csv for the trend
├── R_T_vs_A_H_segmented_<dataset>.csv
├── r_nu_vs_A_H_pivots_<dataset>.csv
├── segmented_curves_<run>.csv                 one row per manifest curve
├── total_average_R_T_<run>.csv
├── leave_one_out_<run>.csv                    the trend refitted with each pooled dataset or experiment left out
├── dataset_diagnostics.csv
├── metadata.toml                              how it was produced: configuration, commit, machine
├── configuration.toml
└── Manifest.toml
plots/<system>/<run>/
├── nu_vs_A.pdf
├── r_nu_vs_A_H.pdf
├── R_T_vs_A_H.pdf
├── r_nu_vs_A_H_segmented_<dataset>.pdf
└── R_T_vs_A_H_segmented_<dataset>.pdf
```

An output follows the same `<quantity>_vs_<abscissa>` rule as an input, so it can be fed back in
without translation, and a table name carries the quantity, the abscissa and the dataset, nothing
more. `<dataset>` is the stem of the input file, `<accession>_<Author>_<year>`, so the archive
accession is in the name of every file of a measurement; the trend's files carry its label,
`systematic_trend`. A label, `Author year`, is the display name and the key that selects a curve
in the manifest, and takes the accession only where two datasets would share it:

| File | Content |
| :--- | :--- |
| `r_nu_vs_A_H_<dataset>.csv` | the ratio extracted point by point from one measurement |
| `R_T_vs_A_H_<dataset>.csv` | the temperature ratio from it |
| `r_nu_vs_A_H_consensus_systematic_trend.csv` | the combined ratio the trend curve was fitted to |
| `r_nu_vs_A_H_segmented_<dataset>.csv` | the fitted ratio, tabulated |
| `R_T_vs_A_H_segmented_<dataset>.csv` | the temperature ratio from the fitted ratio |
| `r_nu_vs_A_H_pivots_<dataset>.csv` | the fit as its joined points |
| `segmented_curves_<run>.csv` | per manifest curve, keyed by its label: accession, segments, pin, span, pairs, coverage, reduced chi-squared and χ² against its expectation, range mean |
| `total_average_R_T_<run>.csv` | ⟨R_T⟩ over each yield distribution reaching the coverage floor, with the yield fraction the curve takes in |
| `leave_one_out_<run>.csv` | the trend refitted with each pooled dataset or experiment left out |
| `dataset_diagnostics.csv` | one row per dataset read, with the qualifiers its retrieval recorded |
| `manifest_<run>.toml` | what the run produced, for a consuming code |
| `metadata.toml` | how it was produced: configuration, commit, machine |

The ratio tables have the headers `A_H,r_nu,r_nu_uncertainty` and `A_H,R_T,R_T_uncertainty`, the
pivot table included; in the point-by-point tables an uncertainty the measurement does not quote is
an empty field, never zero. `segmented_curves_<run>.csv` has the columns `label, accession, kind,
pooled, segments, pinned_at_symmetric_split, first_A_H, last_A_H, pairs, measured_points, coverage,
reduced_chi_squared, chi_squared_over_expectation, weights_imputed, range_mean_R_T,
range_mean_R_T_uncertainty, deviation_autocorrelation`, all dimensionless.
`chi_squared_over_expectation` equals `reduced_chi_squared` for a dataset curve, and the last,
filled for the systematic trend only, is the `ρ` its covariance was formed with.

`total_average_R_T_<run>.csv` has the columns `segmented_curve, mass_yield, R_T, R_T_uncertainty,
R_T_uncertainty_independent_points, R_T_uncertainty_leave_one_out, R_T_leave_one_out_min,
R_T_leave_one_out_max, yield_fraction, R_T_one_more_segment, R_T_two_more_segments`. The three
leave-one-out columns and the last two are filled for the systematic trend only.
`R_T_uncertainty_independent_points` is the approximation of the published tables over the
tabulated points of a curve.

`leave_one_out_<run>.csv` has the columns `dataset_left_out, accession, mass_yield, R_T, segments,
breakpoints, reduced_chi_squared, outcome`, one row per dataset or experiment left out and yield
distribution. The labels of an experiment are joined by ` + `; its accessions, like the
breakpoints, are separated by spaces.

`dataset_diagnostics.csv` has the columns `dataset, accession, points, pairs, first_pair,
last_pair, coverage, outside_physical_range, complement_sum, complement_spread,
without_uncertainties, qualifiers, pooling_weight, flagged, pair_sum_deviation, scale_consistent,
pooled, pooled_with, exclusion_reason, segmented_curve, curve_flagged, curve_flag_reason`.
`outside_physical_range` counts the complete pairs left out because their `r_ν` falls outside
`(0, 1)`. `segmented_curve` reads "segmented curve" or states why there is none: no
complete pair, coverage below the floor, or no fit and the reason.
[`write_results`](@ref) refuses a run whose identifier would make one of those names longer than
the 255 bytes a file system admits.

The manifest, the curve table, the total averages and the leave-one-out refits repeat the
identifier in their names; the per-curve ratio tables do not, and are found through the manifest.
A consuming code stages the whole run directory and selects the manifest by its `manifest_`
prefix, and the token keeps a staged copy attributable after it has left `data/sims/`. A run
directory holds exactly one. [Consuming the output](consuming.md) describes the manifest. The
systematic-trend curve is labelled `systematic_trend`, the same token as its kind.

A run-identifier parameter is the configuration key it came from, by section and key as the file
spells it, `segments.max_segments`. Every key that changes the result of the run is a parameter;
the system is not, since it names the directory the identifier sits in, nor are `[retrieval]`, a
guard on the inputs, and `[output]`, which changes how a number is rendered and not the number. A
key that cannot change this run's result is none of its parameters either: the Gilbert-Cameron
keys of a back-shifted Fermi gas run, `windows_apply_to_datasets` without windows, an input that
is not given. Only the parameters behind the `[domain]` table of the manifest are readable in the
identifier — `ldm`, `avg`, `nZ`, `chg`, `mass`, `dZ0`, and `def` in a Gilbert-Cameron run:

```
avg=charge_resolved_chg=wahl_dZ0=false_h=a72fb8dd_ldm=BSFG_mass=ame2020_nZ=5
```

A path among them, a tabulated charge distribution or mass table, appears as the first eight
hexadecimal digits of the SHA-1 of its spelling relative to the data directory. Every parameter,
the readable ones included, enters `h`: the same eight digits of the SHA-1 of the canonical text of
[`run_parameters`](@ref), one `section.key=value` line per parameter, sorted. A path enters that
text relative to the data directory, so the same inputs staged on another machine give the same
identifier. `metadata.toml` records the parameters under `[identifier.parameters]`, and
[`run_identifier`](@ref) recomputes the identifier from them. The identifier is
to be compared, never parsed: the parameters of a run are read from its `metadata.toml`.
