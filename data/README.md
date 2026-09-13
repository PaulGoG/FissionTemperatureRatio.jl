# Input data

This directory describes the input data; the data itself is held locally and is not part of what
the repository ships. None of it originates with this package, and the MIT licence covering the
source code does not extend to it: each collection carries the terms and the attribution of its
own source, recorded below.

## Layout

| Directory | Content |
|---|---|
| `mass_excess/` | atomic mass evaluation, `Z A symbol D σD`, D and σD in keV |
| `charge_distribution/` | charge polarization and Gaussian dispersion, `A ΔZ rms` |
| `shell_corrections/` | shell corrections for the Gilbert-Cameron systematic, `n S(N) S(Z)` in MeV |
| `multiplicity/<nucleus>/` | experimental prompt neutron multiplicity, `A ν σν` |
| `yield/<nucleus>/` | experimental pre-neutron fragment mass yield, `A Y σY` |

All files are whitespace-separated with a single header line, except the mass excess table, which
has none. Readers take columns by position, not by header text.

The `multiplicity/` and `yield/` directories also hold a `retrieval.toml`, the run record of the
query that produced them, naming every dataset it kept or excluded and why. The readers take only
`.dat` files, so the record sits beside the data without interfering.

## Sources

**`mass_excess/AME2020.ANA`** — the 2020 atomic mass evaluation.
W. J. Huang, M. Wang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030002 (2021);
M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030003 (2021).
Distributed by the Atomic Mass Data Center.

**`charge_distribution/DeltaZ_rms_A.*`** — charge polarization `ΔZ(A)` and the root-mean-square
dispersion `rms(A)` of the isobaric charge distribution, for `U5` (235-U), `PU39` (239-Pu) and
`CF52` (252-Cf).
A. C. Wahl, *Atomic Data and Nuclear Data Tables* **38**, 1–156 (1988).
No table is held for 233-U; see the known defects below.

**`shell_corrections/SZSN.GC`** — shell corrections `S(N)` and `S(Z)`, tabulated against nucleon
number from 11 to 150.
A. Gilbert, A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965), as distributed in
the IAEA Reference Input Parameter Library, segment on level densities, file `Beijing.gc`. Checked
against Table III of that paper: the values and the column order agree.

**`multiplicity/<nucleus>/*.dat`** and **`yield/<nucleus>/*.dat`** — experimental prompt neutron
multiplicity and pre-neutron fragment mass yield against fragment mass, from the EXFOR
Experimental Nuclear Data Library of the IAEA Nuclear Data Services,
<https://www-nds.iaea.org/exfor/>. Each measurement remains the work of its authors and should be
cited as such in any result derived from it; the EXFOR DatasetID that opens each file name
identifies the entry, and the `retrieval.toml` beside it records the reaction code and units.

Retrieved with [ExforFissionData.jl](https://github.com/PaulGoG/ExforFissionData.jl), which writes
this layout directly. To reproduce a directory, clone that package and run its retrieval — from
*its* checkout, not this one — naming this package as the output root:

```
julia --project scripts/retrieve.jl config/U233_nf_yield_A.toml /path/to/FissionTemperatureRatio.jl
```

EXFOR entries are immutable once published, so a configuration and that package reproduce a
retrieval exactly.

## Selecting prompt multiplicity data: the archive coding is not sufficient

Per-fragment multiplicity is nominally marked by the `FRG` tag of the EXFOR reaction code, and a
retrieval that asks for it returns only tagged datasets. For 252-Cf that rule is sound in the
direction it asserts — every tagged set checked here is per fragment — but **unsound as an
exclusion**: of the eight 252-Cf datasets coded `MASS,PR,NU`, without the tag, five are per
fragment and three are per pair, reporting the total multiplicity of the split at 3 to 5 neutrons
per fission.

The discriminator that does work is the data. A per-pair quantity is a property of the split, so
it must be invariant under `A -> A₀ - A`; a per-fragment quantity is a sawtooth whose complementary
values differ by a factor of four or five and sum to about the total multiplicity. Applied to the
untagged 252-Cf sets:

| Dataset | Verdict |
|---|---|
| 14652004 Britt 1964, 23118006 Zeynalov 2011, 23175008 Budtz-Jørgensen 1988, 23268005 Göök 2014, 41689004 Piksaykin 1977 | per fragment, kept |
| 23213012 Mehta 1973, 41720003 Basova 1979, 404200022 Zakharova 1979 | per pair, or the tagged subentry of the same measurement is already held; excluded |

`multiplicity/Cf252_0f/` therefore holds the tagged retrieval plus those five, and carries the
run record of both queries — `retrieval.toml` for the tagged one and `retrieval-untagged.toml` for
the other. This matters beyond tidiness: Budtz-Jørgensen 1988 is the canonical 252-Cf(sf) `ν(A)`
reference and one of the sets the published table averages, so a tag-only selection could not
reproduce it.

The two files that first exposed this carry the header `nuPair errnuPair` from the earlier
retrieval, which named the quantity after the tag rather than after the contents. Readers take
columns by position, so no header ever affected a result.

## Known defects

**No charge distribution table exists for 233-U.** Runs for that nucleus fall back to the average
values `ΔZ = -0.5`, `rms = 0.6`, which the pipeline reports at every run. These are not EXFOR
observables — they are fit parameters of a `Z_p` model — so the retrieval cannot supply them; the
table would have to be transcribed from Wahl. The effect is bounded: at most 3.4 % on `R_T(A_H)`
at the shell minimum with a median of 0.27 %, and far less on the total average, where the charge
treatment largely cancels.

**The charge distribution tables carry no edition record.** `DeltaZ_rms_A.*` were assembled by
hand from two single-column tables and do not say which edition or fit of Wahl they came from.

**Uncertainties are absent from several multiplicity sets**, which is reflected in the
`weights_imputed` field of every fit that used them. Such a set still receives an uncertainty on
its total average, propagated from the yield distribution.
