# Input data

This directory describes the input data; the data itself is held locally and is not part of what
the repository ships. None of it originates with this package, and the MIT licence covering the
source code does not extend to it: each collection carries the terms and the attribution of its
own source, recorded below. Pipeline runs write their output under `data/sims/`, which is not
version-controlled either and is not input.

## Layout

```
data/
└── <system>/                                      Cf252_sf, U233_nth, U235_nth, Pu239_nth
    ├── charge_distribution_vs_A.dat               optional
    ├── nu_vs_A/
    │   ├── retrieval.toml
    │   ├── <accession>_<Author>_<year>.dat
    │   └── subentries/                            the archive subentries read, not data
    ├── Y_vs_A/                                    as nu_vs_A/
    └── TKE_vs_A/                                  as nu_vs_A/
```

The system-independent evaluations are not held here: the atomic mass evaluation, AME2020, and the
Gilbert-Cameron shell corrections, Table III of the 1965 paper, ship with FissionFragmentsDomain.jl,
and a configuration names them `"ame2020"` and `"gilbert_cameron_1965"`. A table of either in the
same layout can be named by its path under `data/` instead. So does Wahl's charge distribution model, which the
shipped configurations use; a tabulated charge distribution is needed only to reproduce the
published extraction exactly.

One directory per system, one subdirectory per measured quantity, one file per measurement: the
directory carries the quantity and the file carries the provenance. A system is named
`<ElementSymbol><A>_<entrance channel>`, and a tabular file `<quantity>_vs_<abscissa>`. An
extension states a format and never a nuclide.

| File | Columns |
|---|---|
| `<system>/charge_distribution_vs_A.dat` | `A dZ sigma_Z` |
| `<system>/nu_vs_A/*.dat` | `A nu nu_uncertainty`, or `A nu` where no uncertainty is quoted |
| `<system>/Y_vs_A/*.dat` | `A Y Y_uncertainty` |
| `<system>/TKE_vs_A/*.dat` | `A TKE TKE_uncertainty`, or `A TKE`; pre-neutron, in MeV |

All files are whitespace-separated with a single header line. **Readers take columns by position, not by
header text**, so a header rename upstream cannot affect a result here — and nothing downstream
may be coupled to header text either.

Each measured-quantity directory also holds a `retrieval.toml`, the run record of the query that
produced it, naming every dataset it kept or excluded and why. The readers take only `.dat` files,
so the record sits beside the data without interfering.

## Sources

**`<system>/charge_distribution_vs_A.dat`** — charge polarization `ΔZ(A)` and the Gaussian
dispersion `σ_Z(A)` of the isobaric charge distribution, held for `U235_nth`, `Pu239_nth` and
`Cf252_sf`: the digitised tables of the published extraction, read only to reproduce it.
A. C. Wahl, *Atomic Data and Nuclear Data Tables* **39**, 1–156 (1988),
[doi:10.1016/0092-640X(88)90016-2](https://doi.org/10.1016/0092-640X(88)90016-2).
These are the per-reaction least-squares fits, not the CYF systematics. The same fits, evaluated
from their parameters rather than digitised, are what FissionFragmentsDomain.jl's Wahl (1988)
model gives.

**`<system>/nu_vs_A/*.dat`**, **`<system>/Y_vs_A/*.dat`** and **`<system>/TKE_vs_A/*.dat`** —
experimental prompt neutron multiplicity, pre-neutron fragment mass yield and pre-neutron mean
total kinetic energy against fragment mass, from the EXFOR
Experimental Nuclear Data Library of the IAEA Nuclear Data Services,
<https://www-nds.iaea.org/exfor/>. Each measurement remains the work of its authors and should be
cited as such in any result derived from it; the EXFOR DatasetID that opens each file name
identifies the entry, and the `retrieval.toml` beside it records the reaction code and units.

[ExforFissionData.jl](https://github.com/PaulGoG/ExforFissionData.jl) is the supported source of
these files and writes this layout directly. To reproduce a directory, clone that package and run
its retrieval from its own checkout, giving the root of this repository as the output root:

```
julia scripts/retrieve.jl config/U233_nth_nu_vs_A.toml /path/to/FissionTemperatureRatio
julia scripts/retrieve.jl config/U233_nth_TKE_vs_A.toml /path/to/FissionTemperatureRatio
```

Check that package out at a named commit before retrieving: the `[run]` table of every
`retrieval.toml` records the revision and version that wrote it, and a run of this package copies
that table, for every input directory it reads, into its `metadata.toml`. The results quoted in the
README rest on the twelve committed configurations of its v0.2.1, `<system>_<observable>.toml`
for the four systems and `nu_vs_A`, `Y_vs_A`, `TKE_vs_A`.

Masses the archive gives as non-integer values, digitised or binned, are interpolated onto the
integers by the retrieval, which records `mass_treatment = "interpolated"` and the count of raw
values. Such a dataset is pooled at the weight of its measured points over its written rows.

A DatasetID has eight digits, or nine where it points within a subentry; the dataset label drops
the leading digits either way. Two files of one author and year keep their DatasetID in
parentheses after the label, so that every label is unique.

EXFOR entries are immutable once published, so a configuration and that package reproduce a
retrieval exactly.

## Selecting prompt multiplicity data: the archive coding is not sufficient

Per-fragment multiplicity is nominally marked by the `FRG` tag of the EXFOR reaction code, but for
252-Cf the tag is sound only in the direction it asserts: several datasets coded `MASS,PR,NU`,
without it, are per fragment, and others are per pair, reporting the total multiplicity of the
split at 3 to 5 neutrons per fission. The discriminator that works is the data. A per-pair
quantity is a property of the split, so it is invariant under `A -> A₀ - A`; a per-fragment
quantity is a sawtooth whose complementary values differ by a factor of four or five and sum to
about the total multiplicity.

The retrieval applies that test and records its verdict on each dataset coded without the tag:

| Dataset | Verdict |
|---|---|
| 23118006 Zeynalov 2011, 23175008 Budtz-Jørgensen 1988, 23268005 Göök 2014 | per fragment, admitted |
| 23213012 Mehta 1973, 41720003 Basova 1979 | per pair; rejected |
| 14652004 Britt 1964, 41689004 Piksaykin 1977 | per fragment by the test, not read: the publications could not be consulted |

Budtz-Jørgensen 1988 is the canonical 252-Cf(sf) `ν(A)` reference and one of the datasets the
published table averages, so a tag-only selection could not reproduce it. Both subentries of
Zakharova 1979, the tagged 404200021 and the untagged 404200022, are refused because their masses
are provisional: the entry names a companion experiment that applied no correction for neutron
emission. For 233-U, 22660006 Nishio 1998 and
41397006 Apalin 1965 are per pair as coded and are rejected in the same way.

## Known defects

**No digitised charge distribution table is held for 233-U.** The published extraction took
`ΔZ = -0.5`, `σ_Z = 0.6` there, which `charge_distribution_file = "mean"` reproduces. The shipped
configuration takes Wahl's 1988 parameters for that reaction instead, from FissionFragmentsDomain.jl.

**The digitised charge distribution tables carry no edition record.** The
`charge_distribution_vs_A.dat` tables were assembled by hand from two single-column tables and do
not say which edition or fit of Wahl they came from. They serve only the exact reproduction of the
published extraction; the shipped configurations do not read them.

**Qualified datasets are used, not corrected.** A dataset whose retrieval record lists `DERIV` or
`SPA` among its reaction-code qualifiers is flagged in the log, in `dataset_diagnostics.csv` and in
the run metadata. Six staged sets carry `SPA`: the 239-Pu `ν(A)` of Basova 1979 and of
Zamyatnin 1979, the 239-Pu `Y(A)` of Walter 1964, and the 235-U `Y(A)` of Romano 2010, Bohn 1969
and Straede 1987, over the last of which two rows of the published table are averaged. Most
thermal-neutron sets carry `MXW`, a Maxwellian-averaged spectrum, which is the entrance channel and
is not flagged.

**Several yield distributions are partial.** Britt 1963 (252-Cf), Bohn 1969 (235-U) and Akimov
1971 (239-Pu) are measured on the light wing alone, Barreau 1985 (23717005) on the far heavy wing
and Dyachenko 1967 (233-U) over twelve heavy masses. With `symmetrize = true` a light-wing yield
stands for its heavy complement; a distribution giving a yield at fewer than
`min_dataset_coverage` of the heavy mass numbers is read but not averaged over, and every total
average states the fraction of its distribution's yield it takes in.

**Uncertainties are absent from several multiplicity datasets.** An unquoted uncertainty is read
as `missing`, never as zero; such points take the median weight of the quoted ones and are counted
in the `weights_imputed` field of every fit that used them and in the `without_uncertainties`
column of the dataset diagnostics. Such a dataset still receives an uncertainty on its total
average, propagated from the fit covariance and the yield distribution.
