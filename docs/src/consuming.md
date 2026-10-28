# Consuming the output

The manifest is the entry point for a code that takes these curves as input, and it is meant to
be read rather than parsed out of file names. It is FissionFragmentsDomain.jl's run record,
written by its `write_temperature_ratio_manifest` and read by its
`read_temperature_ratio_manifest`, and holds nothing else. To stage a run, copy the whole run
directory from `data/sims/<system>/<run identifier>/`; it holds exactly one `manifest_*.toml`,
which `staged_manifest` selects by that prefix. An abridged ²³³U manifest:

```toml
[system]
label = "U233_nth"
notation = "²³³U(nth,f)"
target_A = 233
target_Z = 92
channel = "nth"
reaction = "n,f"
incident_energy_MeV = 2.53e-8
compound_A = 234
compound_Z = 92

[run]
ordinate = "R_T"
abscissa = ["A_H"]
columns = ["A_H", "R_T", "R_T_uncertainty"]

[domain]
level_density_model = "BSFG"
deformed_branch = false
ratio_averaging = "charge_resolved"
excitation_weighted = true
charges_per_mass = 5
charge_model = "Wahl1988(U233T, Z_F = 92, A_F = 234)"
mass_table = "mass_excess_ame2020.dat"
package_version = "0.2.4"
zero_polarization_at_symmetry = false  # absent from manifests written before 0.2.4

[[segmented_curve]]
label = "K. Nishio 1998"
kind = "dataset"                   # or "systematic_trend"
temperature_ratio_file = "R_T_vs_A_H_segmented_22660005_K.Nishio_1998.csv"
multiplicity_ratio_pivots_file = "r_nu_vs_A_H_pivots_22660005_K.Nishio_1998.csv"
accession = "22660005"             # EXFOR dataset identifier; dataset curves only

[[segmented_curve]]
label = "systematic_trend"
kind = "systematic_trend"
temperature_ratio_file = "R_T_vs_A_H_segmented_systematic_trend.csv"
multiplicity_ratio_pivots_file = "r_nu_vs_A_H_pivots_systematic_trend.csv"
```

What else the run reports about each curve is in `segmented_curves_<run identifier>.csv`, one row
per manifest curve keyed by the same label, and in `metadata.toml`.

### The domain must match

`[domain]` records the fragmentation domain the curves were extracted on: the level density model
and its Gilbert-Cameron branch, the averaging and its excitation weighting, the charges per mass,
the charge model, the mass table, whether `ΔZ = 0` was imposed at the symmetric split and the
version of FissionFragmentsDomain.jl. A code that partitions the excitation energy on another
domain does not invert the same relation.
[DeterministicSequentialEmission.jl](https://github.com/PaulGoG/DeterministicSequentialEmission.jl),
which reads these curves, refuses such a curve.

### Read the columns by position, and read `temperature_ratio_file`

`columns` says what the three columns hold; it is not there to be matched against the header
text. `R_T` is not piecewise-linear even where the multiplicity ratio is, because the level
density parameters carry the shell structure into it. It is therefore tabulated at every mass
number, and whatever interpolation a consumer applies reproduces the tabulated value. The pivot
file holds `r_ν`, not `R_T`. The tabulated `R_T_uncertainty` is the square root of the diagonal of
the propagated covariance of the curve. The tabulated points are not independent, so summing
their uncertainties in quadrature understates the uncertainty of any average over them.
`total_average_R_T_<run identifier>.csv` carries the covariance-propagated uncertainty of each
total average first and the independent-points one beside it.

### One curve is taken

The curves are alternatives and are not averaged. Which one describes the system is decided by
running the consuming code with each and comparing its output with experiment. The curve of
`kind = "systematic_trend"`, labelled `systematic_trend`, is meant for use where a dataset is too
sparse or too scattered to resolve the shape.

### The default inverts a per-pair partition exactly

A code that forms `a_L/a_H` for every fragmentation `(A_H, Z_H)` and applies `R_T(A_H)` to each
recovers the measured `r_ν(A_H)` to rounding from a `charge_resolved` curve. From an
excitation-weighted curve it does so when it weights the pairs by their excitation, as the
measured multiplicities do. A `ratio_of_means` curve, the published setting, misses that round
trip by up to 4 × 10⁻³ in `R_T` at the doubly magic heavy fragment.

`test/test_manifest.jl` is a consumer: it reads a run through FissionFragmentsDomain.jl's reader
and nothing else, and checks the contract this page describes.
