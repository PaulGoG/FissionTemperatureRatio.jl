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

A run accepts a retrieval record only at or above `[retrieval] min_package_version`, 0.2.3 in the
shipped configurations, and refuses one that states no `[run] package_version`.

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
README rest on the twelve committed configurations of its v0.2.3, `<system>_<observable>.toml`
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

The retrieval applies that test and records its verdict on each dataset coded without the tag,
with the scale of the reading: the pair sum weighted with the light-fragment yield, against `ν̄`.

| Dataset | Verdict |
|---|---|
| 252-Cf: 23118006 Zeynalov 2011, 23175008 Budtz-Jørgensen 1988, 23268005 Göök 2014, 14652004 Britt 1964, 41689004 Piksaykin 1977 | per fragment, admitted |
| 239-Pu: 22650004 Tsuchiya 2000, 41502006 Batenkov 2004; 235-U: 41502005 Batenkov 2004 | per fragment, admitted |
| 252-Cf: 23213012 Mehta 1973, 41720003 Basova 1979 | per pair; rejected |

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
1971 (239-Pu) are measured on the light wing alone, Barreau 1985 (23717005) on the far heavy wing,
Dyachenko 1967 (233-U) over twelve heavy masses, and Dyachenko 1967 (235-U, 41713004) is a sparse
digitisation of one figure that misses the heavy peak. With `symmetrize = true` a light-wing yield
stands for its heavy complement. Coverage is measured against the primary distribution of the
system, the share of its heavy-fragment yield at the masses a distribution holds; one below
`min_yield_coverage` is read but not averaged over, and every total average states the fraction
of its distribution's yield it takes in.

**Two 252-Cf yield distributions are not held.** The retrieval refuses Vorobiev 2001, 41425015
and 41425016, coded as the inclusive pre-neutron yield: the subentries place them on Fig. 11a
(`NUt = 0`) of Dushin et al., *Nucl. Instrum. Methods A* **516**, 539 (2004),
doi:10.1016/j.nima.2003.09.029, and their mean heavy mass, 146.7 and 146.2 u, lies about 3 u above
the 142.9 to 143.6 u of every inclusive measurement. The `ν(A)` of 41425014, from the same
measurement, has its sawtooth minimum where the others do and is held.

**One 252-Cf multiplicity dataset is kept out of the pool.** Zeynalov 2019 (41739002) is short
of light-fragment neutrons by 15 to 19 %, its pair sum 3.47 against `ν̄` = 3.76, and its `r_ν` lies
above every other set over most of `A_H` = 130 to 160. It is coded per fragment, so the retrieval
forms no pair sum for it and records no scale flag. The shipped configuration excludes it from the
pooled trend by its accession, with the reason; its own curve is still written. The main README
gives the attribution.

**One 235-U multiplicity dataset is kept out of the pool.** Zeynalov 2019 (41738002), from the
same paper and analysis, has its pair sum on the scale of `ν̄` but divides it wrongly between the
wings: over `A_H` = 130 to 150 its light fragments emit 1.02 neutrons against 1.36 to 1.66 in the
comparable sets, its heavy fragments 1.36 against 1.02 to 1.17. The shipped configuration
excludes it from the pooled trend by its accession, with the reason; its own curve is still
written.

**Two 239-Pu multiplicity datasets are off the scale of ν̄.** Tsuchiya 2000 (22650004) and
Batenkov 2004 (41502006) are read per fragment, with pair sums 3.9 ± 0.6 % above and 10.7 ± 1.3 %
below `ν̄`, and their records state `scale_consistent = false`; the run flags them. A uniform scale
cancels in `r_ν`, so both are used. A scale error that varies with mass would not cancel, and the
pair sum cannot see one; by the structural diagnostics of the run, Tsuchiya's pair sums over
`A_H` = 126 to 152 lie within the spread of the other 239-Pu sets, and its larger spread overall
comes from the heavy tail above 154. Tsuchiya's subentry heads its values in percent per fission,
a miscoding the record states as `unit_miscoded`; they are written as tabulated, in neutrons per
fragment. Batenkov's datasets for 239-Pu and 235-U lie on a 4-u grid of odd masses, on which no
mass has its complement at `A₀ = 240` or 236, so neither enters a curve.

**Uncertainties are absent from several multiplicity datasets.** An unquoted uncertainty is read
as `missing`, never as zero; such points take the median weight of the quoted ones and are counted
in the `weights_imputed` field of every fit that used them and in the `without_uncertainties`
column of the dataset diagnostics. Such a dataset still receives an uncertainty on its total
average, propagated from the fit covariance and the yield distribution.
