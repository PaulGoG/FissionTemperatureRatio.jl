# Input data

Every file here is third-party scientific data, redistributed for reproducibility. None of it
originates with this package, and the MIT licence covering the source code does not extend to it:
each collection carries the terms and the attribution of its own source, recorded below.

## Layout

| Directory | Content |
|---|---|
| `mass_excess/` | atomic mass evaluation, `Z A symbol D σD`, D and σD in keV |
| `charge_distribution/` | charge polarization and Gaussian dispersion, `A ΔZ rms` |
| `shell_corrections/` | shell corrections for the Gilbert-Cameron systematic, `n S(N) S(Z)` in MeV |
| `multiplicity/<nucleus>/` | experimental prompt neutron multiplicity, `A ν σν` |

All files are whitespace-separated with a single header line, except the mass excess table, which
has none. Readers take columns by position, not by header text.

## Sources

**`mass_excess/AME2020.ANA`** — the 2020 atomic mass evaluation.
W. J. Huang, M. Wang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030002 (2021);
M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030003 (2021).
Distributed by the Atomic Mass Data Center.

**`charge_distribution/DeltaZ_rms_A.*`** — charge polarization `ΔZ(A)` and the root-mean-square
dispersion `rms(A)` of the isobaric charge distribution, for `U5` (235-U), `PU39` (239-Pu) and
`CF52` (252-Cf).
A. C. Wahl, *Atomic Data and Nuclear Data Tables* **38**, 1–156 (1988).
These files were assembled by hand from two single-column tables and carry no record of which
edition or fit they were taken from. **No table exists for 233-U**, so runs for that nucleus fall
back to the average values `ΔZ = -0.5`, `rms = 0.6`, which the pipeline reports at every run.

**`shell_corrections/SZSN.GC`** — shell corrections `S(N)` and `S(Z)`, tabulated against nucleon
number from 11 to 150.
A. Gilbert, A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965), as distributed in
the IAEA Reference Input Parameter Library, segment on level densities, file `Beijing.gc`.

**`multiplicity/<nucleus>/*.dat`** — experimental prompt neutron multiplicity as a function of
fragment mass, retrieved from the EXFOR Experimental Nuclear Data Library of the IAEA Nuclear Data
Services, <https://www-nds.iaea.org/exfor/>. Each measurement remains the work of its authors and
should be cited as such in any result derived from it.

The first author and year of each measurement are encoded in the file name:

| Fissioning system | Data sets |
|---|---|
| `U233_nf` | Apalin 1965, Nishio 1998, Takamiya 1999 |
| `U235_nf` | Maslin 1967, Boldeman 1971, Nishio 1998, Vorobyev 2010, Al-Adili 2020 |
| `Pu239_nf` | Apalin 1965, Basova 1979, Zamyatnin 1979, Nishio 1995 |
| `Cf252_0f` | Bowman 1963, Nardi 1968, Mehta 1973, Basova 1979, Zakharova 1979, Zamyatnin 1979, Vorobiev 2001, Zeynalov 2011, Göök 2014, Zeynalov 2019, Al-Adili 2020 |

## Known defects

**EXFOR entry identifiers are not recorded.** The files carry the first author and year but not the
EXFOR accession and subentry numbers, so a given file cannot be traced to the exact retrieval it
came from, nor checked against a later revision of that entry. Restoring them is the main reason
to regenerate this directory from a parser rather than to patch it by hand.

**Two files disagree with their archive coding.** `Cf252_0f_nuA_Sh.Zeynalov_2011.dat` and
`Cf252_0f_nuA_A.Goeoek_2014.dat` label their second and third columns `nuPair errnuPair`, and the
EXFOR entries they came from are coded `98-CF-252(0,F)MASS,PR,NU` — per fragment pair — where the
other nine are `MASS,PR/FRG,NU`, per fragment.

Their contents are per fragment. A quantity defined for a fragment pair is a property of the
split and must be invariant under `A -> A0 - A`; these files are not. At `A = 120` they read 3.42
and 3.18 against 0.76 and 0.80 at the complementary mass 132, tracking the uncontested
per-fragment set (Vorobiev, 3.10 and 0.69) closely and failing symmetry by the same factor.
Complementary values sum to 3.87 and 3.96, against a total prompt neutron multiplicity of about
3.76 for the spontaneous fission of 252-Cf; a genuine pair quantity would sum to about 7.5.

Readers take columns by position, so neither the header nor the coding affects any result, and
excluding the two sets shifts the mean temperature ratio by 0.27 %. The discrepancy is recorded
because it matters to whoever re-retrieves this directory: selecting on the `FRG` tag alone would
discard two usable data sets.

**Uncertainties are absent from several sets**, which is reflected in the `weights_imputed` field
of every fit that used them.
