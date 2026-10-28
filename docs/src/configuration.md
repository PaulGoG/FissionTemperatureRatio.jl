# Configuration

A run is driven by one TOML file, `config/<system>.toml`, read by [`load_configuration`](@ref)
against a data directory the caller names. The loader enforces what this page lists: the
sections and keys it knows, the type of every key, its bounds or its enumerated choices, and the
existence of every file named; the rules that tie one key to another, given in the constraint
column, are enforced by the [`Configuration`](@ref) constructor, so a configuration built by hand
is held to them as well. An unknown section or key, a value outside its bounds, a choice it does
not offer, `true` or `false` where a number stands, or a path that does not exist is refused with
a message naming the key, so a run cannot start from a configuration it cannot honour. A key that is absent takes its default; a section marked required
must be present even where every key of it has one. The shipped configurations set every key they
read, so that the file states what the run honoured.

Every key that changes the result is a parameter of the run identifier, under the short key of the
fifth column; [Naming](naming.md) defines the identifier and `metadata.toml` records the
parameters in full.

## `[system]`, required

The keys describe what was irradiated. The reaction code, the fissioning nucleus and the system
label `<symbol><A>_<channel>` are derived from them, and the file is named by that label.

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `target_A` | integer | required | 2 to 400 | — | mass number of the target nuclide |
| `target_Z` | integer | required | 1 to 118, below `target_A` | — | charge number of the target nuclide |
| `channel` | string | required | `"sf"`, `"nth"`, `"nres"`, `"nfast"` | — | entrance channel; fixes the reaction code and the fissioning nucleus |
| `incident_energy` | real | `0.0` | at least 0; 0 for `"sf"` | `E` | incident neutron energy, MeV |

## `[fragmentation]`, required

The heavy-fragment mass range and the charges retained at each mass number. ``A_0`` is the mass
number of the fissioning nucleus; the symmetric split is ``A_0/2``, rounded up for an odd ``A_0``.

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `charges_per_mass` | integer | required | odd, 1 to 15 | `nZ` | charge numbers retained about ``Z_p`` at each mass number |
| `heavy_mass_min` | integer | the symmetric split | from the symmetric split to ``A_0 - 1`` | `AHmin` | smallest heavy-fragment mass number |
| `heavy_mass_max` | integer | required | above `heavy_mass_min`, at most ``A_0 - 1`` | `AHmax` | largest heavy-fragment mass number |
| `charge_distribution_file` | string | `"wahl"` | `"wahl"`, `"mean"`, or the path of a table `A dZ sigma_Z` under the data directory | `chg` | the isobaric charge distribution: Wahl's ``Z_p`` model, the mean values of the published extraction, or a tabulated one |
| `zero_polarization_at_symmetry` | boolean | `false` | not allowed with `"wahl"`, which fixes it | `dZ0` | ``\Delta Z = 0`` at the symmetric split for a table or `"mean"` |

## `[level_density]`, required

The level density parameters and how the inversion averages over the charges.

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `model` | string | required | `"BSFG"`, `"GC"` | `ldm` | back-shifted Fermi gas, or Gilbert-Cameron |
| `mass_excess_file` | string | `"ame2020"` | `"ame2020"`, or the path of a mass-excess table under the data directory | `mass` | the atomic mass evaluation [Huang2021, Wang2021](@cite) |
| `shell_correction_file` | string | `"gilbert_cameron_1965"` | `"gilbert_cameron_1965"`, or the path of a table under the data directory; read by `"GC"` alone | `sc`, a Gilbert-Cameron run only | the shell corrections, Table III of Gilbert and Cameron |
| `deformed_branch` | boolean | `true` | read by `"GC"` alone | `def`, a Gilbert-Cameron run only | eq. (21) of the 1965 paper for the nuclei it counts as deformed; `false` applies eq. (20) throughout |
| `ratio_averaging` | string | `"charge_resolved"` | `"charge_resolved"`, `"ratio_of_means"`, `"mean_of_ratios"` | `avg` | how the per-charge ratio ``a_L/a_H`` enters the inversion |
| `mean_kinetic_energy_file` | string | none | the path of a `TKE_vs_A` dataset under the data directory, listed by the retrieval record beside it; `"charge_resolved"` alone | `TKE` | the pre-neutron ``\langle \mathrm{TKE} \rangle(A)`` in MeV that weights each fragmentation by its mean total excitation |

## `[multiplicity]`, required

The ``\nu(A)`` datasets of the system.

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `subdirectory` | string | required | a directory under the data directory holding `A nu nu_uncertainty` tables | — | where the datasets are read from, `<system>/nu_vs_A` |
| `exclude` | array of tables | none | each `{ accession = "…", reason = "…" }`, the accession an EXFOR dataset identifier of 8 or 9 characters naming a dataset of the directory | `excl` | datasets kept out of the pooled trend; each is still read, fitted where its coverage allows, written and diagnosed |

```toml
[multiplicity]
subdirectory = "Cf252_sf/nu_vs_A"
exclude = [
    { accession = "31845003", reason = "one usable fragment pair in range" },
]
```

An exclusion names its dataset by EXFOR accession. A label, `Author year`, takes the accession of
each dataset as soon as a second one of the same author and year is retrieved, and an exclusion
written against it would then match neither; the accession does not change. An exclusion that
names no dataset of the directory is refused when the configuration is loaded.
`dataset = "<label>"` remains for a tabulation that carries no accession; for a dataset that has
one it is refused, with the accession to write.

## `[yield]`, optional

The pre-neutron mass yield distributions the total average ``\langle R_T \rangle`` is taken over.
Without the section the run reports the mean over the fragment mass range alone.

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `mass_yield_file` | string | required within the section | the path of a `Y_vs_A` dataset under the data directory | `Yf`, without `subdirectory` | the primary distribution: averaged over alone, or the reference the coverage of every distribution of `subdirectory` is measured against |
| `subdirectory` | string | none | a directory under the data directory | `Y`, one table of the directory, the primary and the exclusions | every distribution of the directory to average over in place of the primary alone |
| `exclude` | array of tables | none | as `[multiplicity] exclude`; refused without `subdirectory` | within `Y` | distributions kept out of every average; each is still read and reported |
| `symmetrize` | boolean | `true` | — | `Ysym` | ``Y(A)`` and ``Y(A_0 - A)`` averaged where both are measured, the pre-neutron identity |
| `min_yield_coverage` | real | `0.3` | 0 to 1 | `Ycov` | the share of the primary distribution's heavy-fragment yield a distribution must cover to be averaged over |

## `[segments]`, required

The segmented description of ``r_\nu(A_H)``; every key has a default.

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `max_segments` | integer | `6` | 1 to 12 | `nseg` | largest segment count examined |
| `min_points_per_segment` | integer | `4` | 2 to 50 | `npts` | data points a segment must contain |
| `pin_symmetric_split` | boolean | `true` | needs an even ``A_0`` and `heavy_mass_min` at the symmetric split | `pin` | fix ``r_\nu = 1/2`` at the symmetric split |
| `required_windows` | array of `[first, last]` | none | ordered integer pairs inside `heavy_mass_min:heavy_mass_max` | `win` | windows the trend curve must break in; `[[128, 132]]` in the shipped configurations |
| `windows_apply_to_datasets` | boolean | `false` | — | `wdat` | apply the windows to every dataset curve as well |
| `min_pair_coverage` | real | `0.3` | 0 to 1 | `cov` | the fraction of the range's heavy masses a dataset must pair to offer a curve of its own |
| `autocorrelation_lags` | integer | `4` | 1 to 8 | `lags` | lags of the deviation correlogram the trend's autocorrelation is fitted to |

## `[retrieval]`, optional

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `min_package_version` | string | `"0.2.7"` | a version number | — | the lowest ExforFissionData.jl version whose retrieval records a run accepts; a record below it, or stating none, is refused |

## `[output]`, required

| Key | Type | Default | Constraint | Short key | Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `significant_digits` | integer | `6` | 1 to 15 | — | significant figures in the tabulated output; it changes the rendering and not the number, so it is no parameter of the identifier |
