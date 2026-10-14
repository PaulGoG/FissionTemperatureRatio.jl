# FissionTemperatureRatio.jl

[![CI](https://github.com/PaulGoG/FissionTemperatureRatio.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PaulGoG/FissionTemperatureRatio.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/PaulGoG/FissionTemperatureRatio.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PaulGoG/FissionTemperatureRatio.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://PaulGoG.github.io/FissionTemperatureRatio.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PaulGoG.github.io/FissionTemperatureRatio.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.11%2B-9558B2?logo=julia&logoColor=white)](https://julialang.org)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Code style: JuliaFormatter](https://img.shields.io/badge/code%20style-JuliaFormatter-informational)](https://github.com/domluna/JuliaFormatter.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Extraction of the temperature ratio `R_T = T_L/T_H` of complementary fully accelerated fission
fragments from experimental prompt neutron multiplicity data, described by joined straight
segments. The fragmentation domain, the level density parameters and the relation between `R_T`
and the excitation-energy partition are those of
[FissionFragmentsDomain.jl](https://github.com/PaulGoG/FissionFragmentsDomain.jl), which the
prompt emission codes that read `R_T(A_H)` partition on as well.

The method, its conventions and its equation numbering are those of *Eur. Phys. J. A* **60**, 190
(2024), [doi:10.1140/epja/s10050-024-01375-7](https://doi.org/10.1140/epja/s10050-024-01375-7).
This file is the entry point; the
[documentation](https://PaulGoG.github.io/FissionTemperatureRatio.jl/stable/) sets out the
[method](docs/src/method.md), what has been [validated](docs/src/validation.md) against the paper,
every [configuration](docs/src/configuration.md) key, the contract a
[consuming code](docs/src/consuming.md) reads, and the [naming](docs/src/naming.md) of files and
columns.

```
FissionTemperatureRatio/
├── config/          pipeline configurations, one TOML file per fissioning system
├── data/            input data, held locally, README.md on sources and terms; sims/ for run output
├── scripts/run.jl   pipeline entry point, own environment with the plotting stack
├── src/             configuration, inversion, segmented fit and curve, pipeline, provenance
├── ext/             CairoMakie extension: figures
├── test/            test suite, own environment
├── bench/           benchmarks, own environment
├── docs/            Documenter site, own environment, references.bib; assets.jl renders figures
├── formatter/       JuliaFormatter environment used by check.jl
├── check.jl         format, then test
├── activate.jl      activate and instantiate the root environment
├── CITATION.cff     citation metadata for the package and the method paper
└── Project.toml     dependencies and compatibility bounds
```

The full tree is at the [end of this file](#full-file-tree). Generated output is written to
`data/sims/` and `plots/`, neither of which is version-controlled.

## Environments

The package is not registered: clone the repository and work from its root. Each environment
carries an activation script that activates and instantiates it silently; the first
instantiation resolves and precompiles, and is slow.

```
julia activate.jl
julia scripts/activate.jl
julia test/activate.jl
julia docs/activate.jl
julia bench/activate.jl
```

Every runnable script activates its own environment as its first statement, so every command
below runs as written, with no project flag. For an interactive session on the package,
`julia -i activate.jl`. The formatting environment under `formatter/` is instantiated by
`check.jl` and by CI, not by hand.

Each auxiliary environment takes the package through a relative `[sources]` entry, honoured from
Pkg 1.11, the declared floor, so it runs against the local source and not a registered snapshot;
the activation scripts resolve before instantiating, so that a dependency added to the package
does not leave their manifests stale. FissionFragmentsDomain.jl is not registered either: every
environment takes it from its repository at a tagged release through `[sources]`, with the same
release as the `[compat]` floor, fetched with the git executable (`JULIA_PKG_USE_CLI_GIT`), which
honours the user's git configuration — credentials, `url.insteadOf` rewrites to SSH — where
libgit2 can hang.

Resolved manifests are not version-controlled: a manifest resolved against one Julia version
would break instantiation on the others the package supports. `scripts/run.jl` copies the resolved
manifest of its environment into the run directory as `Manifest.toml`, so a result stays
attributable to the exact dependency versions that produced it.

## Entry points

Run the pipeline for one fissioning nucleus:

```
julia scripts/run.jl config/U233_nth.toml
```

The tables, the manifest and the provenance record go to `data/sims/<system>/<run identifier>/`
and the figures to `plots/<system>/<run identifier>/`; see [the run directory](#the-run-directory).

Run the test suite:

```
julia -e 'include("activate.jl"); using Pkg; Pkg.test()'
```

Apply the formatting gate and then the tests, as CI does:

```
julia check.jl
```

Run the benchmarks:

```
julia bench/benchmarks.jl
```

Build the documentation:

```
julia docs/make.jl
```

Regenerate the figures shown below, from the configuration and the input data:

```
julia docs/assets.jl
```

## Status

| Component | State |
|---|---|
| Mass table, shell corrections, charge distribution, fragmentation domain, level density parameters | FissionFragmentsDomain.jl v0.2.3 |
| Multiplicity ratio | complete |
| Temperature ratio: charge-resolved inversion, excitation-weighted by `⟨TKE⟩(A)` | complete; a `⟨TKE⟩(A)` dataset is configured for each of the four systems |
| Temperature ratio: ratio of means, mean of ratios | complete |
| Retrieval records and reaction-code qualifiers | complete |
| Segmented description with order selection | complete |
| Pipeline, tabulated output, figures, provenance | complete |
| Per-dataset and systematic-trend curves | complete |
| Averaging over a fragment mass yield distribution | complete |
| Reproduction of published total averages | complete for ²³³U, ²⁵²Cf and ²³⁵U, and for ²³⁹Pu over an inferred Y(A); see [Validation](docs/src/validation.md) |
| Run record with the fragmentation domain, for the consuming code | complete |

Versions follow semantic versioning with the pre-1.0 convention: a breaking change, to a
configuration key, to an output file name or to a key of the manifest, raises the minor version,
and a patch release carries none. The [changelog](CHANGELOG.md) marks every breaking change.

## Results at a glance

With the published settings, this package returns the total average temperature ratios of
Table 1 of the paper to within 1.0 % for ²³³U(nth,f), ²⁵²Cf(sf), ²³⁵U(nth,f) and ²³⁹Pu(nth,f),
Apalin's ²³³U measurement aside at 1.9 %, and with the shipped settings, which invert charge by
charge on Wahl's charge model, to within 1.4 %, Apalin aside at 2.1 %. For ²³⁹Pu the yield
distribution is inferred by agreement, the table's caption naming none. The paper did not use
this package: its jointed segments were drawn by hand through the multiplicity ratio of each
dataset. Here the number of segments and their breakpoints are selected by a criterion, and the
trend that combines several datasets is this package's own construction, which is why the
agreement is at the per-cent level. [Validation](docs/src/validation.md) states exactly what has
been checked against the paper and what has not.

Experimental input — prompt neutron multiplicity `ν(A)`, fragment mass yields `Y(A)` and the mean
total kinetic energy `⟨TKE⟩(A)` — comes from the IAEA EXFOR archive through
[ExforFissionData.jl](https://github.com/PaulGoG/ExforFissionData.jl), the supported route, which
retrieves fission observables and writes them as tabulated files with a record of every dataset it
kept or excluded. This package reads those files; it does not query the archive itself.

![Temperature ratio](docs/src/assets/temperature_ratio.png)

The temperature ratio of ²⁵²Cf(sf): every measurement held, in grey, and the systematic trend
through the fourteen that are pooled, with its uncertainty. Above the symmetric split the light
fragment is the hotter of the pair, and the ratio crosses unity near the most probable
fragmentation.

![Segment selection](docs/src/assets/segment_selection.gif)

One measurement from each of the four configured systems, fitted with an increasing number of
joined segments. Points carry the uncertainty propagated from the multiplicity data where the
archive quotes one. The breakpoints are not placed by hand: for each order they are searched
exhaustively over the abscissae the data occupies, and the order itself is chosen by the Bayesian
information criterion, which prices every added segment and every added breakpoint. The four
systems do not stop at the same order — five, six, four and five — which is why they are shown
together.

The systematic trend of each shipped configuration, with `⟨R_T⟩` over the primary yield
distribution, the uncertainty propagated from the fit covariance and, beside it, the jackknife
over the refits with each pooled dataset or experiment left out:

| System | `R_T` ± `R_T_uncertainty` | `R_T_uncertainty_leave_one_out` |
|---|---|---|
| ²⁵²Cf(sf) | 1.0821 ± 0.0159 | 0.0085 |
| ²³⁵U(nth,f) | 1.1230 ± 0.0153 | 0.0222 |
| ²³⁹Pu(nth,f) | 1.0388 ± 0.0269 | 0.0257 |
| ²³³U(nth,f) | 1.0724 ± 0.0936 | 0.1381 |

These rest on inputs retrieved by ExforFissionData.jl v0.2.7; v0.2.8 writes the same tables and
differs in its records alone.

## Figures

The plotting stack is a weak dependency, so `using FissionTemperatureRatio` loads in about a
second and pulls in no graphics. Figures come from an extension that loads with CairoMakie:

```julia
using FissionTemperatureRatio, CairoMakie
```

`scripts/run.jl` runs in the environment under `scripts/`, which carries CairoMakie, so a
pipeline run always writes its figures. A caller using the package as a library gets the tables,
the manifest and `run_metadata` without it, and figures through `write_figures` with it.

## Input data

The input data is not shipped with the package. It is third-party scientific data — experimental
prompt neutron multiplicity, fragment mass yield and kinetic energy measurements — held locally
under the terms of its own sources. [`data/README.md`](data/README.md) records what each file is
and where it comes from, the known defects of the staged sets, and the attribution of every
exclusion the shipped configurations make. The atomic mass evaluation, the Gilbert-Cameron shell
corrections and Wahl's charge distribution parameters are FissionFragmentsDomain.jl's, which
ships them.

Place it under `data/` in the layout that file describes:

```
data/
└── <system>/                           Cf252_sf, U233_nth, U235_nth, Pu239_nth
    ├── charge_distribution_vs_A.dat    A dZ sigma_Z         (optional)
    ├── nu_vs_A/                        A nu nu_uncertainty, or A nu
    ├── Y_vs_A/                         A Y Y_uncertainty    (optional)
    └── TKE_vs_A/                       A TKE TKE_uncertainty, TKE in MeV, pre-neutron (optional)
```

One directory per system, one subdirectory per measured quantity, one file per measurement, each
with one header line; readers take columns by position, not by header text. The measurements come
from [ExforFissionData.jl](https://github.com/PaulGoG/ExforFissionData.jl), the supported source of
`ν(A)`, `Y(A)` and `⟨TKE⟩(A)`, which writes exactly this layout, every file name carrying its
EXFOR DatasetID, beside a `retrieval.toml` run record naming every dataset the query kept or
excluded. Run its retrieval from its own checkout, with the root of this repository as the second
argument, the output root:

```
julia scripts/retrieve.jl config/U233_nth_nu_vs_A.toml /path/to/FissionTemperatureRatio
```

A run accepts a retrieval record only at or above `[retrieval] min_package_version`, `"0.2.7"` by
default, the first release whose records name the datasets of one experiment, and refuses a
record below it or stating no version. The reaction-code qualifiers a record carries (`DERIV`,
`SPA`, `preliminary`) and a scale it finds off `ν̄` are reported and flagged in the log, in
`dataset_diagnostics.csv` and in the run metadata; no correction and no exclusion follows from a
flag.

Without the data the pipeline cannot run, and the tests that need measurements are reported as
skipped in the test summary. The rest of the suite — the inversion on the shipped fragmentation
domains, the configuration, a synthetic run and its run record — runs regardless, on the mass
table and shell corrections FissionFragmentsDomain.jl ships.

## Configuration

Every run is driven by a TOML file under `config/`, which fixes the fissioning nucleus, the
fragmentation range, the level density model, the datasets and their exclusions, the yield
distributions, the segment search and the precision of the tabulated output. The parser enforces
the types, enumerated choices and bounds its comments document, refuses unknown sections and
keys, and fails naming the offending key, so a run cannot start from a configuration it cannot
honour. Every key, with its type, bounds, default and run-identifier token, is listed under
[Configuration](docs/src/configuration.md).

A system is declared by what was irradiated, not by what fissions: the target, the entrance
channel and the incident energy. The reaction code, the fissioning nucleus and the system label
`<symbol><A>_<channel>` — `Cf252_sf`, `U235_nth` — are all derived from them, so a label cannot
contradict the nuclide it names, and two incident energies of the same target are two systems
rather than one. The two sections that fix the physics, from the shipped ²³³U configuration:

```toml
[fragmentation]
charges_per_mass = 5                   # odd, 1 to 15
heavy_mass_min = 117                   # at least A0/2
heavy_mass_max = 159                   # below A0
charge_distribution_file = "wahl"      # "wahl" | "mean", or a table A dZ sigma_Z under data/
zero_polarization_at_symmetry = false  # dZ = 0 at A0/2 for a table or "mean"

[level_density]
model = "BSFG"                         # "BSFG" | "GC"
mass_excess_file = "ame2020"           # the evaluation FissionFragmentsDomain.jl ships, or a table
shell_correction_file = "gilbert_cameron_1965"  # Table III as FissionFragmentsDomain.jl ships it, or a table
deformed_branch = true                 # "GC": eqs. (20) and (21); false for eq. (20) throughout
ratio_averaging = "charge_resolved"    # "charge_resolved" | "ratio_of_means" | "mean_of_ratios"
mean_kinetic_energy_file = "U233_nth/TKE_vs_A/21981008_P.Geltenbort_1985.dat"
```

The fragmentation domain is built by FissionFragmentsDomain.jl from these keys; `"wahl"` resolves
Wahl's `Zₚ` model for the system. `ratio_averaging = "charge_resolved"`, the default, solves for
`R_T` in the relation between `R_T` and `E*_H/TXE` reduced over the charge distribution
fragmentation by fragmentation — the exact inverse of a prompt emission code that partitions every
pair with its own `a_L/a_H` — and `mean_kinetic_energy_file` weights every fragmentation by its
mean total excitation `Q + E*_CN − ⟨TKE⟩(A_H)`. `"ratio_of_means"` is the closed form of the
published extraction, which the published settings reproduce exactly over the digitised charge
tables; [the method](docs/src/method.md) gives why the default changed. The shipped configurations
take `⟨TKE⟩(A)` and `Y(A)` from one primary experiment per system — Göök 2014 for ²⁵²Cf,
Al-Adili 2016 for ²³⁵U, Wagemans 1984 for ²³⁹Pu and Geltenbort 1985 for ²³³U — and the run
metadata names the accessions used.

To add a fissioning system, place its `ν(A)` datasets in `data/<system>/nu_vs_A/` and copy a
configuration to `config/<system>.toml`.

## Several measurements of one system

Each dataset is carried through the whole chain on its own: its own `r_ν`, its own segmented fit,
its own `R_T`, its own total average. A run offers one curve per measurement that reaches the
coverage floor and one more, the systematic trend, fitted to the combination of those pooled. They
are alternatives offered to a prompt emission code, not an ensemble to be averaged: a consuming
code takes exactly one, and which one describes reality is settled downstream, by comparing the
multiplicity distributions and yields that code produces against experiment.

The datasets of one system disagree far beyond their quoted uncertainties, by ten or twenty times
the median quoted uncertainty for ²⁵²Cf, and that disagreement shapes the rest.

- **Pooling.** The trend is fitted to a curve combined mass number by mass number with an
  additional between-dataset variance, so that the weights are nearly equal and the uncertainty of
  the combination reflects the disagreement instead of hiding it. Datasets of one experiment,
  which their retrieval records name, enter the pool as one.
- **Uncertainty.** The error of a pooled dataset is taken as correlated along the mass axis, with
  a coefficient estimated from the pooled datasets themselves, and the covariance of the trend is
  the sandwich form for such errors; taken as independent, the uncertainty would be low by a
  factor of two to four. What the correlation does not hold, a constant offset of one dataset, is
  measured by refitting the trend with each dataset or experiment left out, and the jackknife over
  the refits is written beside the propagated uncertainty.
- **Admission.** No dataset is rejected for being far from the others, since none is consistent
  with any other in units of its quoted uncertainty. Every dataset is described by structural
  diagnostics in `dataset_diagnostics.csv`, and one is kept out of the pool only when the
  configuration names it by accession and says why, or when its record marks it as superseded by
  a republication. The shipped configurations exclude Zeynalov 2019 for ²⁵²Cf and for ²³⁵U,
  Nishio 1995 for ²³⁹Pu, and Batenkov 2004 for ²³⁵U and ²³⁹Pu; an excluded dataset is still read,
  fitted, written and diagnosed.

[The method](docs/src/method.md) gives the construction and the numbers of the four shipped
systems, and [`data/README.md`](data/README.md) the attribution of each exclusion.

## The run directory

`run_pipeline(configuration)` returns an `ExtractionResult` and writes nothing;
`write_results(result, directory)` writes the tables and the manifest into `directory` and refuses
one that already holds files. `scripts/run.jl` names that directory
`data/sims/<system>/<run identifier>/`, moves an existing run of the same identifier aside as
`<run identifier>#1`, `#2`, …, writes the figures to `plots/<system>/<run identifier>/` through
`write_figures`, and adds the provenance record `metadata.toml` — the library's `run_metadata`,
the configuration and script paths, a timestamp, the git commit and the platform — beside copies
of the configuration, as `configuration.toml`, and of the resolved manifest of the scripts
environment, as `Manifest.toml`.

| File | Content |
|---|---|
| `manifest_<run identifier>.toml` | the run record a consuming code reads; exactly one per directory |
| `r_nu_vs_A_H_<dataset>.csv`, `R_T_vs_A_H_<dataset>.csv` | the ratios extracted point by point from one measurement |
| `r_nu_vs_A_H_consensus_systematic_trend.csv` | the combined ratio the trend curve was fitted to |
| `r_nu_vs_A_H_segmented_<dataset>.csv`, `R_T_vs_A_H_segmented_<dataset>.csv` | the fitted ratio at every mass number and the temperature ratio from it; `…_systematic_trend.csv` for the trend |
| `r_nu_vs_A_H_pivots_<dataset>.csv` | the fit as its joined points; `…_systematic_trend.csv` for the trend |
| `segmented_curves_<run identifier>.csv` | one row per manifest curve, keyed by its label: segments, pin, span, pairs, coverage, reduced chi-squared, range mean |
| `total_average_R_T_<run identifier>.csv` | `⟨R_T⟩` of every segmented curve over every yield distribution |
| `leave_one_out_<run identifier>.csv` | the systematic trend refitted with each pooled dataset or experiment left out; written where the run made the refits |
| `dataset_diagnostics.csv` | one row per dataset read, with its reaction-code qualifiers |
| `metadata.toml`, `configuration.toml`, `Manifest.toml` | the provenance record, written by the script |

`<dataset>` is the stem of the input file, `<accession>_<Author>_<year>` as the retrieval writes
it, so `R_T_vs_A_H_segmented_22660005_K.Nishio_1998.csv` traces to its EXFOR entry and to its
input by name alone. The label of a curve stays `Author year`, with the accession added only where
two datasets would share it: it is the display name and the key a consuming code selects a curve
by. The run identifier is DrWatson's `savename` over the parameters behind the manifest's
`[domain]` table, with `dZ0`, and `h`, the first eight hexadecimal digits of the SHA-1 of the
canonical text of every result-changing configuration key; `metadata.toml` records the parameters
in full under `[identifier.parameters]` and the text `h` was taken over. For ²³³U:

```
avg=charge_resolved_chg=wahl_dZ0=false_h=a72fb8dd_ldm=BSFG_mass=ame2020_nZ=5
```

[Naming](docs/src/naming.md) lists the columns of every table, and
[Configuration](docs/src/configuration.md) the short key of every configuration key.

## Consuming the output

The manifest is the entry point for a code that takes these curves as input, and it is meant to
be read rather than parsed out of file names. It is FissionFragmentsDomain.jl's run record,
written by its `write_temperature_ratio_manifest` and read by its
`read_temperature_ratio_manifest`: the system, the fragmentation domain the curves were extracted
on and, per segmented curve, its label, its kind, its two files and, for a dataset curve, its
accession. A consuming code stages the whole run directory and selects the manifest by its
prefix. Four things to know before writing against it:

- **The domain must match.** A code that partitions the excitation energy on another domain does
  not invert the same relation, and [DeterministicSequentialEmission.jl](https://github.com/PaulGoG/DeterministicSequentialEmission.jl), which reads
  these curves, refuses such a curve.
- **Read the columns by position, and read `temperature_ratio_file`.** `R_T` is not
  piecewise-linear where `r_ν` is, so it is tabulated at every mass number; the pivot file holds
  `r_ν`. The tabulated uncertainties are not independent, so summing them in quadrature
  understates the uncertainty of any average over them.
- **The curves are alternatives, not an ensemble.** Take one; `systematic_trend` where a dataset
  is too sparse or too scattered to resolve the shape.
- **The default inverts a per-pair partition exactly.** A `ratio_of_means` curve, the published
  setting, misses that round trip by up to 4 × 10⁻³ in `R_T` at the doubly magic heavy fragment.

[Consuming the output](docs/src/consuming.md) gives an abridged manifest and the contract in
full; `test/test_manifest.jl` is a consumer that checks it.

## Method

For each fragment pair the prompt neutron multiplicity ratio is identified with the excitation
energy ratio, and the fragment level densities are taken in the Fermi-gas regime, so that one
fragmentation with `ρ = a_L/a_H` gives the heavy fragment `E*_H/TXE = 1/(1 + ρ R_T²)`. The
measured `r_ν = ν_H/(ν_L + ν_H)` is resolved by mass alone, and `R_T(A_H)` is the root of

```
r_ν(A_H) = Σ_Z p(Z, A_H) ⟨TXE⟩_Z / (1 + ρ_Z R_T²)  /  Σ_Z p(Z, A_H) ⟨TXE⟩_Z,
```

with `⟨TXE⟩_Z = Q(A_H, Z) + E*_CN − ⟨TKE⟩(A_H)`, which involves no fit and no prompt emission
calculation. With one effective ratio `R_a` in place of `ρ_Z` it reduces to the closed form
`R_T = [(1 - r_ν)/(R_a r_ν)]^(1/2)` of the published method. Because the ratio extracted point by
point is scattered, it is `r_ν` that is described by a continuous piecewise-linear function whose
segment count and breakpoints are selected from the data by the Bayesian information criterion,
and `R_T` follows by the inversion above. Constraints that are exact are enforced rather than
fitted: `r_ν` is pinned to one half at the symmetric split, for every curve whose data reach it,
and the fitted curve may not leave `(0, 1)`, outside which the relation is undefined.

[The method](docs/src/method.md) sets this out in full, with references, including the
combination of datasets, the uncertainty of the trend and why the charge-resolved inversion
replaced the ratio of means as the default.

## How to cite

FissionTemperatureRatio.jl is written by Paul-Adrian Gogîță. For the method, cite
*Eur. Phys. J. A* **60**, 190 (2024),
[https://doi.org/10.1140/epja/s10050-024-01375-7](https://doi.org/10.1140/epja/s10050-024-01375-7);
for the implementation, cite this repository by its URL,
<https://github.com/PaulGoG/FissionTemperatureRatio.jl>. `CITATION.cff` at the root carries both
in machine-readable form, the paper as the preferred citation. As BibTeX:

```bibtex
@software{FissionTemperatureRatio,
  author  = {Gogîță, Paul-Adrian},
  title   = {FissionTemperatureRatio.jl},
  version = {0.3.1},
  year    = {2026},
  url     = {https://github.com/PaulGoG/FissionTemperatureRatio.jl}
}
```

## Licensing

The source code is under the MIT licence in `LICENSE`.

That licence does not extend to the contents of `data/`, none of which originates with this
package. The atomic mass evaluation, the charge distribution systematics, the shell corrections
and every prompt neutron multiplicity and fragment mass yield measurement are third-party
scientific data, held locally and not redistributed here — `data/` carries only its own README.
That file records what each input is, where it came from and how it should be cited; any result
derived from a measurement should cite that measurement.

## Full file tree

<details>
<summary>Every tracked directory and entry file</summary>

```
FissionTemperatureRatio/
├── .github/                        continuous integration and dependency updates
│   ├── dependabot.yml
│   └── workflows/                  CI.yml (tests, coverage, docs) and format.yml (formatting gate)
├── .JuliaFormatter.toml            formatting rules, enforced by check.jl and CI
├── activate.jl                     activate and instantiate the root environment
├── CHANGELOG.md                    notable changes
├── CITATION.cff                    citation metadata for the package and the method paper
├── check.jl                        pre-commit: format, then test
├── LICENSE                         MIT, covering the source code only
├── Project.toml                    dependencies and compatibility bounds
├── bench/                          benchmark suite, own environment
│   ├── activate.jl
│   ├── benchmarks.jl
│   └── Project.toml
├── config/                         pipeline configurations, one per fissioning system
│   ├── Cf252_sf.toml
│   ├── Pu239_nth.toml
│   ├── U233_nth.toml
│   └── U235_nth.toml
├── data/                           input data, held locally and not version-controlled
│   ├── README.md                   what each input file is, where it comes from, its terms
│   └── sims/                       run output, data/sims/<system>/<run identifier>/
├── docs/                           Documenter site, own environment
│   ├── activate.jl
│   ├── assets.jl                   regenerates the figures the README and the documentation show
│   ├── make.jl
│   ├── Project.toml
│   └── src/                        index, method, validation, configuration, consuming, naming, reference, references; assets/ for figures
│       └── references.bib          DOI-keyed bibliography of the method and validation pages
├── ext/                            weak-dependency extensions
│   └── FissionTemperatureRatioCairoMakieExt.jl   publication figures, loaded with CairoMakie
├── formatter/                      JuliaFormatter environment, version-bounded
│   ├── activate.jl
│   └── Project.toml
├── scripts/                        pipeline entry point, own environment
│   ├── activate.jl
│   ├── Project.toml
│   └── run.jl                      pipeline entry point
├── src/
│   ├── autocorrelation.jl          errors of pooled datasets along the mass axis: correlogram, decay
│   ├── configuration.jl            TOML configuration, parsed and validated; model builders
│   ├── consensus.jl                combining datasets, and the diagnostics describing them
│   ├── extracted_curve.jl          one extracted curve: fit, R_T, covariance; its averages
│   ├── FissionTemperatureRatio.jl  module definition and public interface
│   ├── kinetic_energy.jl           ⟨TKE⟩(A) onto the heavy-mass range; retrieval records
│   ├── multiplicity_ratio.jl       ν(A) input and the multiplicity ratio
│   ├── pipeline.jl                 the run, from input to tabulated output and run record
│   ├── plotting.jl                 figure interface; the implementation is in ext/
│   ├── provenance.jl               run identification and metadata
│   ├── segmented_fit.jl            continuous piecewise-linear regression
│   ├── temperature_ratio.jl        the inversion over a tabulated r_ν
│   └── yields.jl                   Y(A) directories and the independent-points total average
└── test/                           test suite, own environment
    ├── activate.jl
    ├── Project.toml
    ├── runtests.jl
    ├── fixtures.jl                 locates the staged input data; the tests that need it are skipped without it
    ├── test_<module>.jl            one file per source module
    ├── test_experiments.jl         datasets of one experiment pooled as one
    ├── test_manifest.jl            a consumer: reads a run through FissionFragmentsDomain.jl alone
    ├── test_published_values.jl    the total averages against Tables 1 and 2 of the paper
    └── test_quality.jl             Aqua, JET, ExplicitImports
```

</details>
