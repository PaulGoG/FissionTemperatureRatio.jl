```@meta
CurrentModule = FissionTemperatureRatio
```

# Naming

A quantity has exactly one name, and that name is used everywhere it appears: Julia identifier,
struct field, configuration key, directory name, file name, column header, figure label, prose.
The spelling changes register between contexts — `heavy_mass_max` in a configuration, `A_H` in a
path — but the *quantity* does not.

The corollary is the rule that decides the hard cases: **a name says what the thing is, not what
was convenient to type.** `R_T_vs_A_H` is a statement; `R_T_A_H` is a list of symbols. `Data` as a
suffix carries nothing. A column called `value` forces the reader to know the file name in order
to know what the file holds.

This page is the vocabulary as this package applies it. The retrieval that supplies its
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
read from, while [`FragmentationDomain`](@ref) has a field `A_H_range`, mirroring the equations.

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
| ν(*A*) | `nu` | `ν` | prompt neutron multiplicity per fragment |
| *Y*(*A*) | `Y` | `Y` | pre-neutron fragment mass yield |
| *r*_ν | `r_nu` | `r_ν` | multiplicity ratio, ν_H/(ν_L + ν_H) |
| *R*_a | `R_a` | `R_a` | level density parameter ratio, a_L/a_H |
| *R*_T | `R_T` | `R_T` | temperature ratio, T_L/T_H |

`rms` is retired. It named a root-mean-square of nothing stated; the quantity is the Gaussian
dispersion of the isobaric charge distribution, and the field writes σ_Z.

### Uncertainty

One spelling, one position: **`<quantity>_uncertainty`**, immediately after the quantity it
belongs to. `A nu nu_uncertainty`, never `A nu errnu`, and never a bare `uncertainty` column that
leaves the reader to work out which quantity it belongs to. In Julia the symbol is used instead —
`σν`, `σY` — because there the quantity it belongs to is not in doubt.

## Julia identifiers

| Kind | Case | Example |
| :--- | :--- | :--- |
| module, type, abstract type | `CamelCase` | [`SegmentedCurve`](@ref), [`LevelDensityModel`](@ref) |
| function, macro | `lowercase_snake_case` | [`level_density_parameter`](@ref), [`system_label`](@ref) |
| variable, field, keyword | `lowercase_snake_case` | `charges_per_mass`, `significant_digits` |
| constant | `SCREAMING_SNAKE_CASE` | [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref) |

Unicode is used where an identifier *is* a symbol of the governing equations — `ΔZ`, `σ_Z`, `ν`,
`Zₚ`, `A₀`, `R_T` — and only there, so that the implementation reads like the formulation. The
exported API stays ASCII-typeable: every public function can be called, and every keyword argument
of one given, without a compose key. That is why [`plot_multiplicities`](@ref) takes `A_0` while
the field it is passed is `A₀`.

Verb prefixes are a small closed set — `read_`, `write_`, `build_`, `load_`, `run_`, `is_` — and a
function that fits none of them is named for what it returns. `read_<thing>` must return the type
named `<Thing>`, which is the rule that decides [`read_mass_excess_table`](@ref) against
`read_mass_excess`: the type is [`MassExcessTable`](@ref), so the reader carries `table` too.
`load_<thing>` reads *and* validates, which is why the configuration entry point is
[`load_configuration`](@ref) and not `read_configuration`.

Type names carry no `Data` suffix: a type holding multiplicities is [`Multiplicity`](@ref), one
holding mass yields is [`MassYield`](@ref), one holding the charge polarization and dispersion is
[`ChargeDistribution`](@ref). `Table` survives only where the thing genuinely is a lookup table
keyed by nucleon number. No adjective stands in for a noun, and no abstract noun stands in for a
description: the piecewise-linear description of one measurement is a [`SegmentedCurve`](@ref),
not a `Parameterization` — "segmented" says how, "parameterized" says nothing a reader can act on.

The result type of a package is `<Verb>Result`, **one per package**. Here it is
[`ExtractionResult`](@ref), named for what this package does: it extracts a temperature ratio.
`PipelineResult` was not a verb, and `RunResult` would collide, on `using` both, with the same
name in a code that consumes these curves.

Abbreviations are permitted only where the field itself uses them — `TKE`, `TXE`, `BSFG`, `GC`,
`EXFOR`, `AME`, `RIPL`, `sf`, `BIC` — and never as `err`, `param`, `calc`, `val`, `idx`, `num`,
`cfg`, `pts`. The one place short forms are unavoidable is a run identifier, and they are
collected there in a single table; see below.

## Configuration files

A configuration is named `<system>.toml`, the token that names the data directories it reads and
the output directory it writes. A system carrying more than one configuration would name them
`<system>_<case>.toml`; none does. The test suite asserts the name against the system each file
declares, so a file cannot drift from what it runs.

Sections are `[system]`, `[fragmentation]`, `[level_density]`, `[multiplicity]`, `[yield]`,
`[segments]` and `[output]`, and **a key is never prefixed with its section's name**:
`[level_density] model`, not `level_density_model`.

A key's comment states three things and nothing else: what the key is, in one short line; the
enumerated choices where it is constrained; and the bounds and the unit. No physics narrative, no
run rationale — those belong in the documentation. The loader enforces exactly what the comment
declares, so a run cannot start from a configuration it cannot honour.

**A key that can only be redundant or wrong is not in the file.** The reaction code follows from
the entrance channel through [`CHANNEL_REACTION`](@ref), so it is not written down; the fissioning
nucleus follows from the target and the channel, so it is not written down either.

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
`U235_nth`, `U235_nres`, `U233_nth`, `Pu239_nth`. `0f` is a reaction code, not a name, and it
belongs in the run record where the reaction code already is. The channel is also what
distinguishes a thermal from a resonance run of one system, which `n,f` alone cannot.

[`system_label`](@ref) returns that token. The typeset form for a figure or a caption —
`²⁵²Cf(sf)`, `²³³U(nth,f)` — is [`system_notation`](@ref), a separate function: one name may not
mean both.

## Directories, files and headers

```
data/
├── reference/                        system-independent evaluations
│   ├── mass_excess_ame2020.dat
│   └── shell_corrections_gilbert_cameron.dat
└── <system>/                         Cf252_sf, U233_nth, U235_nth, Pu239_nth
    ├── charge_distribution_vs_A.dat
    ├── nu_vs_A/
    │   ├── retrieval.toml            the run record of the retrieval
    │   └── <accession>_<Author>_<year>.dat
    └── Y_vs_A/
        └── <accession>_<Author>_<year>.dat
```

One directory per system, one subdirectory per measured quantity, one file per measurement. The
directory carries the quantity and the file carries the provenance, which is what makes the
archive accession number recoverable from a figure legend.

A tabular file is named `<quantity>_vs_<abscissa>[_<abscissa>]`. `vs` is what makes the name a
statement rather than a list of symbols. Extensions state a format and nothing else: `.dat`
whitespace-separated, `.csv` comma-separated, `.toml` metadata. **A file extension never carries a
nuclide**, because the nuclide is already in the path.

One header line, ASCII, abscissae first, then the ordinate, then its uncertainty:

```
A nu nu_uncertainty
A Y Y_uncertainty
A dZ sigma_Z
A_H,R_T,R_T_uncertainty
```

**No column is named `value`.** Naming the column for its quantity is what frees the reader from
having to know the file name to know what the file is holding, and it matters most for the files
that cross a package boundary.

Readers take columns **by position**, not by header text. That is what makes a header rename a
no-op for code, upstream and downstream alike, and it must stay that way: every reader in this
package says so in its docstring, and nothing may be changed to a lookup by name. The one
evaluation held as distributed, the atomic mass table, has no header line at all — which is the
same rule seen from the other side.

## Output files and run identifiers

A run is one directory, and the directory carries the run identifier. Tables go under the system
in `data/sims/`, figures under the system in `plots/`:

```
data/sims/<system>/<run>/
├── manifest_<run>.toml                        what the run produced, for a consuming code
├── r_nu_vs_A_H_<dataset>.csv
├── R_T_vs_A_H_<dataset>.csv
├── r_nu_vs_A_H_consensus_systematic_trend.csv
├── r_nu_vs_A_H_segmented_<label>.csv
├── R_T_vs_A_H_segmented_<label>.csv
├── r_nu_vs_A_H_pivots_<label>.csv
├── total_average_R_T.csv
├── dataset_diagnostics.csv
├── metadata.toml                              how it was produced: configuration, commit, machine
├── configuration.toml
└── Manifest.toml
plots/<system>/<run>/
├── nu_vs_A.pdf
├── r_nu_vs_A_H.pdf
├── R_T_vs_A_H.pdf
├── r_nu_vs_A_H_segmented_<label>.pdf
└── R_T_vs_A_H_segmented_<label>.pdf
```

An output follows the same `<quantity>_vs_<abscissa>` rule as an input, so it can be fed back in
without translation, and a table name carries the quantity, the abscissa and the label of the
curve or dataset, nothing more:

| File | Content |
| :--- | :--- |
| `r_nu_vs_A_H_<dataset>.csv` | the ratio extracted point by point from one measurement |
| `R_T_vs_A_H_<dataset>.csv` | the temperature ratio from it |
| `r_nu_vs_A_H_consensus_systematic_trend.csv` | the combined ratio the trend curve was fitted to |
| `r_nu_vs_A_H_segmented_<label>.csv` | the fitted ratio, tabulated |
| `R_T_vs_A_H_segmented_<label>.csv` | the temperature ratio from the fitted ratio |
| `r_nu_vs_A_H_pivots_<label>.csv` | the fit as its joined points |
| `total_average_R_T.csv` | ⟨R_T⟩ over each yield distribution |
| `dataset_diagnostics.csv` | one row per dataset read |
| `manifest_<run>.toml` | what the run produced, for a consuming code |
| `metadata.toml` | how it was produced: configuration, commit, machine |

The manifest is the one file inside the directory whose name repeats the identifier. A consuming
code stages the whole run directory and selects the manifest by its `manifest_` prefix, and the
token keeps a staged copy attributable after it has left `data/sims/`. A run directory holds
exactly one.

**A run-identifier token is the configuration key it came from.** Spelling every key out in full
would put some two hundred characters into every directory name, so the keys are abbreviated —
and abbreviated in one place, [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref), rather than at the point of
use. That table is keyed by the key's dotted path in the configuration file,
`segments.min_segment_span => minspan`, so an entry names exactly one key.
[`run_identifier`](@ref) refuses a key with no entry rather than inventing a token, and the test
suite asserts that every key it uses has one.

Every key that changes the result is a token. The system is not: it names the directory the
identifier sits in. `significant_digits` is not either, since it changes how a number is rendered
and not the number. A value that is a list or a path — the required windows, the exclusion list,
the charge distribution file, the yield directory — does not reduce to a token and enters as the
first eight hexadecimal digits of the SHA-1 of its canonical spelling, paths taken relative to the
data directory so that the same inputs staged on another machine give the same identifier. The
value itself is written in full into `metadata.toml` under `[identifier.hashed]`.
