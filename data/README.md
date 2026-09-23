# Input data

This directory describes the input data; the data itself is held locally and is not part of what
the repository ships. None of it originates with this package, and the MIT licence covering the
source code does not extend to it: each collection carries the terms and the attribution of its
own source, recorded below. Pipeline runs write their output under `data/sims/`, which is not
version-controlled either and is not input.

## Layout

```
data/
├── reference/                                     the system-independent evaluations
│   ├── mass_excess_ame2020.dat
│   └── shell_corrections_gilbert_cameron.dat
└── <system>/                                      Cf252_sf, U233_nth, U235_nth, Pu239_nth
    ├── charge_distribution_vs_A.dat
    ├── nu_vs_A/
    │   ├── retrieval.toml
    │   └── <accession>_<Author>_<year>.dat
    └── Y_vs_A/
        ├── retrieval.toml
        └── <accession>_<Author>_<year>.dat
```

One directory per system, one subdirectory per measured quantity, one file per measurement: the
directory carries the quantity and the file carries the provenance. A system is named
`<ElementSymbol><A>_<entrance channel>`, and a tabular file `<quantity>_vs_<abscissa>`. An
extension states a format and never a nuclide.

| File | Columns |
|---|---|
| `reference/mass_excess_ame2020.dat` | `Z A symbol mass_excess mass_excess_uncertainty`, the last two in keV |
| `reference/shell_corrections_gilbert_cameron.dat` | `n S_N S_Z` in MeV |
| `<system>/charge_distribution_vs_A.dat` | `A dZ sigma_Z` |
| `<system>/nu_vs_A/*.dat` | `A nu nu_uncertainty`, or `A nu` where no uncertainty is quoted |
| `<system>/Y_vs_A/*.dat` | `A Y Y_uncertainty` |

All files are whitespace-separated with a single header line, except the mass excess table, which
is held as the evaluation distributes it and has none. **Readers take columns by position, not by
header text**, so a header rename upstream cannot affect a result here — and nothing downstream
may be coupled to header text either.

Each measured-quantity directory also holds a `retrieval.toml`, the run record of the query that
produced it, naming every dataset it kept or excluded and why. The readers take only `.dat` files,
so the record sits beside the data without interfering.

## Sources

**`reference/mass_excess_ame2020.dat`** — the 2020 atomic mass evaluation.
W. J. Huang, M. Wang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030002 (2021);
M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030003 (2021).
Distributed by the Atomic Mass Data Center.

**`<system>/charge_distribution_vs_A.dat`** — charge polarization `ΔZ(A)` and the Gaussian
dispersion `σ_Z(A)` of the isobaric charge distribution, held for `U235_nth`, `Pu239_nth` and
`Cf252_sf`.
A. C. Wahl, *Atomic Data and Nuclear Data Tables* **39**, 1–156 (1988),
[doi:10.1016/0092-640X(88)90016-2](https://doi.org/10.1016/0092-640X(88)90016-2).
These are the per-reaction least-squares fits, not the CYF systematics.

**`reference/shell_corrections_gilbert_cameron.dat`** — shell corrections `S_N` and `S_Z`,
tabulated against nucleon number from 11 to 150.
A. Gilbert, A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965), as distributed in
the IAEA Reference Input Parameter Library, segment on level densities, file `Beijing.gc`. Checked
against Table III of that paper: the values and the column order agree.

**`<system>/nu_vs_A/*.dat`** and **`<system>/Y_vs_A/*.dat`** — experimental prompt neutron
multiplicity and pre-neutron fragment mass yield against fragment mass, from the EXFOR
Experimental Nuclear Data Library of the IAEA Nuclear Data Services,
<https://www-nds.iaea.org/exfor/>. Each measurement remains the work of its authors and should be
cited as such in any result derived from it; the EXFOR DatasetID that opens each file name
identifies the entry, and the `retrieval.toml` beside it records the reaction code and units.

[ExforFissionData.jl](https://github.com/PaulGoG/ExforFissionData.jl) is the supported source of
these files and writes this layout directly. To reproduce a directory, clone that package and run
its retrieval from its own checkout, giving the root of this repository as the output root:

```
julia scripts/retrieve.jl config/U233_nth_nu_vs_A.toml /path/to/FissionTemperatureRatio
```

A DatasetID has eight digits, or nine where it points within a subentry; the dataset label drops
the leading digits either way. Two files of one author and year keep their DatasetID in
parentheses after the label, so that every label is unique.

EXFOR entries are immutable once published, so a configuration and that package reproduce a
retrieval exactly.

## Selecting prompt multiplicity data: the archive coding is not sufficient

Per-fragment multiplicity is nominally marked by the `FRG` tag of the EXFOR reaction code, and a
retrieval that asks for it returns only tagged datasets. For 252-Cf that rule is sound in the
direction it asserts — every tagged dataset checked here is per fragment — but **unsound as an
exclusion**: of the eight 252-Cf datasets coded `MASS,PR,NU`, without the tag, five are per
fragment and three are per pair, reporting the total multiplicity of the split at 3 to 5 neutrons
per fission.

The discriminator that does work is the data. A per-pair quantity is a property of the split, so
it must be invariant under `A -> A₀ - A`; a per-fragment quantity is a sawtooth whose complementary
values differ by a factor of four or five and sum to about the total multiplicity. Applied to the
untagged 252-Cf datasets:

| Dataset | Verdict |
|---|---|
| 14652004 Britt 1964, 23118006 Zeynalov 2011, 23175008 Budtz-Jørgensen 1988, 23268005 Göök 2014, 41689004 Piksaykin 1977 | per fragment, kept |
| 23213012 Mehta 1973, 41720003 Basova 1979, 404200022 Zakharova 1979 | per pair, or the tagged subentry of the same measurement is already held; excluded |

`Cf252_sf/nu_vs_A/` therefore holds the tagged retrieval plus those five, and carries the run
record of both queries — `retrieval.toml` for the tagged one and `retrieval_untagged.toml` for the
other. The second record names the observable it *queried*, `nu_bar_vs_A`, while sitting in the
`nu_vs_A` directory, because that is what happened: the query asked for the per-pair quantity and
the five datasets it contributed were verified per fragment before being kept.

This matters beyond tidiness: Budtz-Jørgensen 1988 is the canonical 252-Cf(sf) `ν(A)` reference
and one of the datasets the published table averages, so a tag-only selection could not reproduce
it.

## Known defects

**No evaluated charge distribution table is held for 233-U.** An evaluated
`charge_distribution_vs_A.dat` of the same construction as the other three can be staged under
`data/U233_nth/`. Until it is, runs for that nucleus fall back to `ΔZ = -0.5`, `σ_Z = 0.6`, the
yield-weighted means of the evaluated tables themselves, and the pipeline reports the fallback at
every run. These are not EXFOR observables — they are fit parameters of a `Z_p` model — so the
retrieval cannot supply them. The effect is bounded at 3.4 % on `R_T(A_H)` at the shell minimum,
with a median of 0.27 %.

**The charge distribution tables carry no edition record.** The `charge_distribution_vs_A.dat`
tables were assembled by hand from two single-column tables and do not say which edition or fit of
Wahl they came from.

**Uncertainties are absent from several multiplicity datasets.** An unquoted uncertainty is read
as `missing`, never as zero; such points take the median weight of the quoted ones and are counted
in the `weights_imputed` field of every fit that used them and in the `without_uncertainties`
column of the dataset diagnostics. Such a dataset still receives an uncertainty on its total
average, propagated from the fit covariance and the yield distribution.
