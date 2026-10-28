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

The system-independent evaluations are not held here. The atomic mass evaluation, AME2020, the
Gilbert-Cameron shell corrections, Table III of the 1965 paper, and Wahl's charge distribution
model, which the shipped configurations use, ship with FissionFragmentsDomain.jl. A configuration
names the first two `"ame2020"` and `"gilbert_cameron_1965"`, or names a table of either in the
same layout by its path under `data/`. A tabulated charge distribution is needed only to
reproduce the published extraction exactly.

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

All files are whitespace-separated with a single header line. Readers take columns by position,
not by header text, so a header rename upstream cannot affect a result here. Nothing downstream
may be coupled to header text either.

Each measured-quantity directory also holds a `retrieval.toml`, the run record of the query that
produced it, naming every dataset it kept or excluded and why. The readers take only `.dat` files,
so the record sits beside the data without interfering.

A run accepts a retrieval record only at or above `[retrieval] min_package_version`, 0.2.7 in the
shipped configurations and by default. That is the first version whose records name the datasets
of one experiment. A run refuses a record that states no `[run] package_version`. A lower floor
can still be set explicitly, and a run whose input directories were written by different versions
says so in the log. Below 0.2.7 the datasets of an experiment would be pooled side by side rather
than as one, giving for ²⁵²Cf `⟨R_T⟩` = 1.0878 in place of 1.0821. Inputs retrieved before 0.2.3
can moreover carry the ²³⁹Pu multiplicities of Tsuchiya 2000 (22650004) off the scale of their
subentry, and the ²⁵²Cf yields 41425015 and 41425016, which that release refuses. A directory
holding no retrieval record, such as a tabulation of one's own, is not judged. A `⟨TKE⟩(A)`
dataset is accepted only if a run record beside it lists it. Its record — the parser revision that
wrote it, the entry, a SHA-1 of the whole record — goes into the run metadata.

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
total kinetic energy against fragment mass, from the EXFOR Experimental Nuclear Data Library of
the IAEA Nuclear Data Services, <https://www-nds.iaea.org/exfor/>. Each measurement remains the
work of its authors and should be cited as such in any result derived from it. The EXFOR
DatasetID that opens each file name identifies the entry, and the `retrieval.toml` beside it
records the reaction code and units.

[ExforFissionData.jl](https://github.com/PaulGoG/ExforFissionData.jl) is the supported source of
these files and writes this layout directly. To reproduce a directory, clone that package and run
its retrieval from its own checkout, giving the root of this repository as the output root:

```
julia scripts/retrieve.jl config/U233_nth_nu_vs_A.toml /path/to/FissionTemperatureRatio
julia scripts/retrieve.jl config/U233_nth_TKE_vs_A.toml /path/to/FissionTemperatureRatio
```

Check that package out at a named commit before retrieving. The `[run]` table of every
`retrieval.toml` records the revision and version that wrote it, and a run of this package copies
that table, for every input directory it reads, into its `metadata.toml`. The results quoted in the
README rest on the twelve committed configurations of its v0.5.0, `<system>_<observable>.toml`
for the four systems and `nu_vs_A`, `Y_vs_A`, `TKE_vs_A`.

Every table that retrieval writes is byte-identical to that of every release from v0.2.4 on; the
records differ. Those of v0.2.7 name the datasets of one experiment and carry the `preliminary`
qualifier. Those of v0.2.8 add the qualifier `multiplicity_from_masses` to Nishio 1995 (23012008)
and restate the curation of Tsuchiya 2000 (22650004). Those of v0.4.0 give every yield the keys
`unit_subentry` and `yield_sum`, restate its `normalisation` in the unit of the file, add the
isomer keys to every entry and reword the conventions `abscissa_resolution` and
`duplicate_abscissa`. Those of v0.5.0 list, per directory, the datasets with lines left out as
compilation defects under `archive_defect_warning` (14101003 in `Cf252_sf/TKE_vs_A`, which no
configuration reads). Of the keys that changed after v0.2.7, this package reads none; it reads
`correlated_with`, `correlation_relation` and the qualifiers, which v0.2.7 introduced.

Masses the archive gives as non-integer values, digitised or binned, are interpolated onto the
integers by the retrieval, which records `mass_treatment = "interpolated"` and the count of raw
values. Such a dataset is pooled at the weight of its measured points over its written rows.

A DatasetID has eight digits, or nine where it points within a subentry; the dataset label drops
the leading digits either way. Two files of one author and year keep their DatasetID in
parentheses after the label, so that every label is unique.

**Datasets of one experiment.** A retrieval record names, under `correlated_with`, the other
datasets of the experiment a dataset belongs to, and under `correlation_relation` how they are
related. Pooled datasets that name one another enter the systematic trend as one curve, formed
according to the stated relation; the rule for each relation is given on the method page of the
documentation, under Combining datasets. Each member still offers its own segmented curve. The
staged records give Basova 1979 and Zamyatnin 1979 as alternative analyses of one experiment, for
²⁵²Cf (41720002 and 41694002) and for ²³⁹Pu (41720004 and 41694003). Their tables differ by 0.34
and 0.36 neutrons rms. No staged `ν(A)` set is superseded.

The retrieval stores beside each table the archive subentries it read and records the parser
revision, so a retrieval can be reproduced from the stored subentries even if the archive revises
an entry.

## Selecting prompt multiplicity data: the archive coding is not sufficient

Per-fragment multiplicity is nominally marked by the `FRG` tag of the EXFOR reaction code, but for
²⁵²Cf the tag is sound only in the direction it asserts. Several datasets coded `MASS,PR,NU`,
without it, are per fragment, and others are per pair, reporting the total multiplicity of the
split at 3 to 5 neutrons per fission. The discriminator that works is the data. A per-pair
quantity is a property of the split, so it is invariant under `A -> A₀ - A`; a per-fragment
quantity is a sawtooth whose complementary values differ by a factor of four or five and sum to
about the total multiplicity.

The retrieval applies that test and records its verdict on each dataset coded without the tag,
with the scale of the reading: the pair sum weighted with the light-fragment yield, against `ν̄`.

| Dataset | Verdict |
|---|---|
| ²⁵²Cf: 23118006 Zeynalov 2011, 23175008 Budtz-Jørgensen 1988, 23268005 Göök 2014, 14652004 Britt 1964, 41689004 Piksaykin 1977 | per fragment, admitted |
| ²³⁹Pu: 22650004 Tsuchiya 2000, 41502006 Batenkov 2004; ²³⁵U: 41502005 Batenkov 2004 | per fragment, admitted |
| ²⁵²Cf: 23213012 Mehta 1973, 41720003 Basova 1979 | per pair; rejected |

Budtz-Jørgensen 1988 is the canonical ²⁵²Cf(sf) `ν(A)` reference and one of the datasets the
published table averages, so a tag-only selection could not reproduce it. Both subentries of
Zakharova 1979, the tagged 404200021 and the untagged 404200022, are refused because their masses
are provisional: the entry names a companion experiment that applied no correction for neutron
emission. For ²³³U, 22660006 Nishio 1998 and 41397006 Apalin 1965 are per pair as coded and are
rejected in the same way.

## Known defects

**No digitised charge distribution table is held for ²³³U.** The published extraction took
`ΔZ = -0.5`, `σ_Z = 0.6` there, which `charge_distribution_file = "mean"` reproduces. The shipped
configuration takes Wahl's 1988 parameters for that reaction instead, from
FissionFragmentsDomain.jl. The published-settings run of the test suite takes the mean values for
²³³U whether or not a table is staged beside its inputs. Over a table, Apalin's row of Table 1
comes back 2.3 % off in place of 1.9 %.

**The digitised charge distribution tables carry no edition record.** The
`charge_distribution_vs_A.dat` tables were assembled by hand from two single-column tables and do
not say which edition or fit of Wahl they came from. They serve only the exact reproduction of the
published extraction; the shipped configurations do not read them.

**Qualified datasets are used, not corrected.** A dataset whose retrieval record lists `DERIV`,
`SPA` or `preliminary` among its qualifiers is flagged in the log, in `dataset_diagnostics.csv` and
in the run metadata; no exclusion follows from the flag. Six staged sets carry `SPA`: the ²³⁹Pu
`ν(A)` of Basova 1979 and of Zamyatnin 1979, the ²³⁹Pu `Y(A)` of Walter 1964, and the ²³⁵U `Y(A)`
of Romano 2010, Bohn 1969 and Straede 1987. Two rows of the published table are averaged over the
last of these. Most thermal-neutron sets carry `MXW`, a Maxwellian-averaged spectrum, which is the
entrance channel and is not flagged. `preliminary` marks the archive's status code `PRELM`, or a
publication that calls its results preliminary. Five staged sets carry it: the `ν(A)` of Batenkov
2004 for ²³⁵U (41502005) and ²³⁹Pu (41502006), both excluded already and forming no pair, the ²³⁵U
`ν(A)` of Vorobyev 2010 (41516012), which is pooled, and the ²³⁵U `Y(A)` 417380041 and 417380042.
For a multiplicity read by the retrieval's complement test the record also states the scale of
the reading, the yield-weighted pair sum against `ν̄`. A dataset it finds off that scale is used,
since a uniform scale cancels in `r_ν`. It is flagged in the log, in the `pair_sum_deviation` and
`scale_consistent` columns and in the run metadata. A scale error that varies with mass does not
cancel, and the pair sum cannot see one; the structural diagnostics of the dataset bear on it.

**Several yield distributions are partial.** Britt 1963 (²⁵²Cf), Bohn 1969 (²³⁵U) and Akimov
1971 (²³⁹Pu) are measured on the light wing alone, Barreau 1985 (23717005) on the far heavy wing,
and Dyachenko 1967 (²³³U) over twelve heavy masses. Dyachenko 1967 (²³⁵U, 41713004) is a sparse
digitisation of one figure that misses the heavy peak. With `symmetrize = true` a light-wing yield
stands for its heavy complement. Coverage is measured against the primary distribution of the
system, as the share of its heavy-fragment yield at the masses a distribution holds. A
distribution below `min_yield_coverage` is read but not averaged over, and every total average
states the fraction of its distribution's yield it takes in.

**Two ²⁵²Cf yield distributions are not held.** The retrieval refuses Vorobiev 2001, 41425015
and 41425016, coded as the inclusive pre-neutron yield: the subentries place them on Fig. 11a
(`NUt = 0`) of Dushin et al., *Nucl. Instrum. Methods A* **516**, 539 (2004),
doi:10.1016/j.nima.2003.09.029, and their mean heavy mass, 146.7 and 146.2 u, lies about 3 u above
the 142.9 to 143.6 u of every inclusive measurement. The `ν(A)` of 41425014, from the same
measurement, has its sawtooth minimum where the others do and is held.

The shipped configurations exclude the three datasets below from the pooled trend, each by
accession and with its reason; each is still fitted and written.

**Zeynalov 2019 (41739002) is kept out of the ²⁵²Cf pool.** Its `r_ν` lies 0.03 to 0.08 above
the consensus of the other fourteen sets over `A_H` = 130 to 160, where their spread is 0.02, and
above every one of them at 25 of those 31 masses. The departure is a deficit of light-fragment
neutrons. Weighted with Göök's yields, `⟨ν_L⟩` is 1.71 against 2.01 to 2.12 in Göök 2014,
Budtz-Jørgensen 1988, Zeynalov 2011 and Al-Adili 2020, and `⟨ν_H⟩` is 1.76 against 1.65 to 1.77.
The pair sum is 3.47, 8 % below `ν̄`, and the heavy fragments would emit half the neutrons. A
scale error confined to one wing does not cancel in `r_ν`, and the retrieval forms no pair sum for
a set coded per fragment, so nothing flags it. The subentry gives pre-neutron masses from a twin
ionisation chamber with one neutron detector, corrected for energy loss, pulse-height defect and
background, with statistical uncertainties only. It states no selection in `TKE` or angle and no
correction for the detection geometry. The paper, *EPJ Web Conf.* **211**, 04003 (2019),
[doi:10.1051/epjconf/201921104003](https://doi.org/10.1051/epjconf/201921104003), presents the
result as a check of revised analysis software on 0.5 × 10⁶ coincidences of data taken earlier at
IRMM. The software corrects the fragment kinetic energies for neutron emission. Those energies fix
the fragment velocities that carry the neutron spectrum into the fragment frame, and a deficit
that grows with fragment velocity is what an error there produces. Nothing shows it to be the
better measurement, and the measurement published from IRMM in 2011 (23118006) is pooled. With it
in the pool the ²⁵²Cf trend takes 4 segments, with `χ²/dof` = 2.16 and `⟨R_T⟩` = 1.0756; without
it, 4 segments, 2.02 and 1.0821.

**Zeynalov 2019 (41738002) is kept out of the ²³⁵U pool.** It is the ²³⁵U result of the same
paper and analysis. Here the total is right and its division between the wings is not. Weighted
with the Al-Adili 2016 yields over `A_H` = 130 to 150, `⟨ν_L⟩` and `⟨ν_H⟩` are 1.02 and 1.36, against 1.36 to 1.66 and 1.02 to 1.17 in Fraser 1966, Nishio 1998,
Al-Adili 2020 and Vorobyev 2010, so `r_ν` is 0.571 against 0.413 to 0.436 (Apalin 1965: 0.486),
while the pair sum, 2.39, is on the scale of `ν̄`. At single masses, `ν(112)` = 3.43 against 1.75
to 2.31 in Nishio, Al-Adili and Vorobyev, and `ν(124)` = 1.24 against 0.36 and 0.47 in the last
two. Over that range its `r_ν` lies 0.05 to 0.16 above the consensus of the other seven sets,
whose spread is 0.04, and above every one of them at 19 of the 21 masses. The subentry gives
pre-neutron masses and statistical uncertainties only. With it in the pool the ²³⁵U trend takes 6
segments, with `χ²/dof` = 0.19 and `⟨R_T⟩` = 1.0743; without it, 5 segments, 0.44 and 1.1230.

**Nishio 1995 (23012008) is kept out of the ²³⁹Pu pool.** Its `ν(m*)` is not a neutron
measurement. The paper, *J. Nucl. Sci. Technol.* **32**, 404 (1995),
[doi:10.1080/18811248.1995.9731725](https://doi.org/10.1080/18811248.1995.9731725), measures the
velocities and the energies of both fragments. It takes the pre-neutron mass `m*` from the two
velocities and the post-neutron mass `m` from the energy and velocity of a fragment, and gives
`ν(m*) = m* − m`. No neutron is detected, whereas the other pooled ²³⁹Pu sets count neutrons. It
states its own total, `⟨ν_total⟩` = 3.2 ± 0.1, "about 10 % larger than the evaluated value of
JENDL-3", and sets the difference aside as not affecting its argument; the EXFOR table quotes no
uncertainties. Nor is it independent of a set already pooled: the correction for the plasma delay
of the detectors is fixed through momentum conservation with the `ν(m*)` of Apalin et al. The
excess lies on the heavy wing alone. Weighted with the Wagemans 1984 yields over `A_H` = 130 to
150, `⟨ν_H⟩` is 1.68 against 1.20 to 1.48 in Fraser 1966, Tsuchiya 2000, Apalin 1965 and Zamyatnin
1979, at a `⟨ν_L⟩` of 1.53 against 1.45 to 1.57; the pair sum is 3.21, 12 % above `ν̄` = 2.88. Mass
by mass its `ν_H` lies within the range of the other four at `A_H` = 126 to 135 and above all four
at every mass from 136 to 154. An excess on one wing does not cancel in `r_ν`. The rows of the
published table that rest on it, as a `ν(A)` dataset and through its yield distribution, are
unaffected. With it in the pool the ²³⁹Pu trend takes 6 segments, with `χ²/dof` = 0.54 and
`⟨R_T⟩` = 1.0113; without it, 4 segments, 0.73 and 1.0388.

**Two ²³⁹Pu multiplicity datasets are off the scale of ν̄.** Tsuchiya 2000 (22650004) and
Batenkov 2004 (41502006) are read per fragment, with pair sums 3.9 ± 0.6 % above and 10.7 ± 1.3 %
below `ν̄`, and their records state `scale_consistent = false`; the run flags them. A uniform scale
cancels in `r_ν`, so both are used. A scale error that varies with mass would not cancel, and the
pair sum cannot see one; by the structural diagnostics of the run, Tsuchiya's pair sums over
`A_H` = 126 to 152 lie within the spread of the other ²³⁹Pu sets, and its larger spread overall
comes from the heavy tail above 154. Tsuchiya's subentry heads its values in percent per fission,
a miscoding the record states as `unit_miscoded`; they are written as tabulated, in neutrons per
fragment. Batenkov's datasets for ²³⁹Pu and ²³⁵U lie on a 4-u grid of odd masses, on which no
mass has its complement at `A₀ = 240` or 236, so neither enters a curve. Each gives 21 masses. The
paper, *AIP Conf. Proc.* **769**, 1003 (2005),
[doi:10.1063/1.1945175](https://doi.org/10.1063/1.1945175), presents them as preliminary, with
statistical uncertainties only, and reports substantial deviations from earlier ²³⁵U data in the
symmetric region and at the edges of the distribution. They are not paired against interpolated
complements, and the shipped configurations exclude both by accession, with that reason; the
exclusions change no result.

**Uncertainties are absent from several multiplicity datasets.** An unquoted uncertainty — an
absent field, or one that is zero, negative or `NaN` — is read as `missing`, never as zero; such
points take the median weight of the quoted ones and are counted in the `weights_imputed` field
of every fit that used them and in the `without_uncertainties` column of the dataset diagnostics.
Such a dataset still receives an uncertainty on its total average, propagated from the fit
covariance and the yield distribution. A field that is not numeric is refused with the file and
the line, so that one malformed uncertainty does not pass for a dataset that quotes none.
