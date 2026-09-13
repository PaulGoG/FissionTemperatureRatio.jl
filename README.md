# FissionTemperatureRatio.jl

Extraction of the temperature ratio `R_T = T_L/T_H` of complementary fully accelerated fission
fragments from experimental prompt neutron multiplicity data, and its parameterization by joined
straight segments.

```
FissionTemperatureRatio/
├── activate.jl                     activate and instantiate the root environment
├── CHANGELOG.md                    notable changes, and what was corrected in the rewrite
├── Project.toml                    dependencies and compatibility bounds
├── bench/                          benchmark suite, own environment
│   ├── activate.jl
│   ├── benchmarks.jl
│   └── Project.toml
├── config/                         pipeline configurations, one per fissioning nucleus
│   ├── Cf252_0f.toml
│   ├── Pu239_nf.toml
│   ├── U233_nf.toml
│   └── U235_nf.toml
├── data/                           input data, held locally and not version-controlled
│   └── README.md                   what each input file is, where it comes from, its terms
├── docs/                           Documenter site, own environment
│   ├── activate.jl
│   ├── make.jl
│   ├── Project.toml
│   └── src/
├── scripts/
│   └── run.jl                      pipeline entry point
├── src/
│   ├── configuration.jl            TOML configuration, parsed and validated
│   ├── FissionTemperatureRatio.jl  module definition and public interface
│   ├── fragmentation.jl            fragmentation range and isobaric charge distribution
│   ├── level_density.jl            level density parameter systematics
│   ├── mass_data.jl                mass excess input
│   ├── multiplicity_ratio.jl       ν(A) input and the multiplicity ratio
│   ├── pipeline.jl                 the run, from input to tabulated output
│   ├── plotting.jl                 publication figures
│   ├── provenance.jl               run identification and metadata
│   ├── segmented_fit.jl            continuous piecewise-linear regression
│   └── temperature_ratio.jl        level density parameter ratio and R_T
└── test/                           test suite, own environment
```

Generated output is written to `results/` and `plots/`, neither of which is version-controlled.

## Environments

Each environment carries an activation script that activates and instantiates it silently. The
first instantiation resolves and precompiles, and is slow.

Resolved manifests are not version-controlled: the package supports a range of Julia versions and
a manifest is resolved against one of them, so committing one would break instantiation on the
others. The dependency versions a run actually used are recorded in its metadata, so a result
stays attributable to the code that produced it.

```
julia --project -e 'include("activate.jl")'
julia --project=test -e 'include("test/activate.jl")'
julia --project=docs -e 'include("docs/activate.jl")'
julia --project=bench -e 'include("bench/activate.jl")'
```

## Entry points

Run the pipeline for one fissioning nucleus:

```
julia --project scripts/run.jl config/U233_nf.toml
```

Run the test suite:

```
julia --project -e 'using Pkg; Pkg.test()'
```

Run the benchmarks:

```
julia --project=bench bench/benchmarks.jl
```

Build the documentation:

```
julia --project=docs docs/make.jl
```

## Input data

The input data is not shipped with the package. It is third-party scientific data — an atomic mass
evaluation, charge distribution systematics, shell corrections, and experimental prompt neutron
multiplicity measurements — held locally for development and testing under the terms of its own
sources. `data/README.md` records what each file is and where it comes from.

Place it under `data/` in the layout that file describes:

```
data/
├── charge_distribution/   A ΔZ rms
├── mass_excess/           Z A symbol D σD
├── multiplicity/<case>/   A ν σν
└── shell_corrections/     n S(N) S(Z)
```

Without it the pipeline cannot run, and the tests that exercise the systematics on real nuclides
are skipped with a warning; the rest of the suite uses synthetic inputs and runs regardless.

## Configuration

Every run is driven by a TOML file under `config/`, which fixes the fissioning nucleus, the
fragmentation range, the level density prescription, the segment search and the output layout.
The parser enforces the types, enumerated choices and bounds its comments document, and fails
naming the offending key, so a run cannot start from a configuration it cannot honour.

To add a fissioning nucleus, place its ν(A) data sets in a directory under `data/multiplicity/`
as whitespace-separated `A ν σν` tables with one header line, and copy a configuration.

Two level density prescriptions are available. The back-shifted Fermi gas is the default; setting
`prescription = "GC"` selects Gilbert-Cameron, which for fission fragments returns markedly larger
parameters away from closed shells. Running both bounds a systematic uncertainty that the
propagated experimental uncertainties do not cover.

## Method

For each fragment pair the prompt neutron multiplicity ratio is identified with the excitation
energy ratio, and the fragment level densities are taken in the Fermi-gas regime. With
`r_ν = ν_H/(ν_L + ν_H)` and the level density parameter ratio `R_a = a_L/a_H`,

```
R_T = [(1 - r_ν) / (R_a r_ν)]^(1/2),
```

which involves no fit and no prompt emission calculation. Because the ratio extracted point by
point is scattered, it is `r_ν` that is parameterized — by a continuous piecewise-linear function
whose segment count and breakpoints are selected from the data by the Bayesian information
criterion — and `R_T` follows by the relation above.

Constraints that are exact are enforced rather than fitted: the charge polarization vanishes at
the symmetric split, where the two fragments are the same nuclide; `r_ν` is pinned to one half
there; and the parameterization may not leave `(0, 1)`, outside which the relation above is
undefined. Together these make `R_T(A₀/2) = 1` hold exactly, as an outcome rather than an
imposition.

A run produces one parameterization per experimental data set, plus a systematic-trend curve
fitted through all of them with the minimum at the heavy magic fragment placed rather than fitted.
These are alternatives for a prompt emission code to choose between, not an ensemble to be
averaged: the data sets of one fissioning nucleus can disagree well beyond their quoted
uncertainties, and which curve describes reality is settled downstream, by comparing the
multiplicity distributions and yields the code produces against experiment.

`docs/src/method.md` sets this out in full, with references.

## Validation status

What has been checked, and what has not. Stated explicitly because the method is published and a
reader's first question is whether this reproduces it.

Verified:

- The exact identities at the symmetric split. `R_a = 1` and `R_T = 1` hold to machine precision
  where the two fragments are the same nuclide, as an outcome of the construction rather than an
  imposition, and a run reports it if they do not.
- The algebra of the extraction. `R_T` inverts the excitation energy partition exactly, the
  uncertainty propagation matches its closed form, and the parameterization is continuous at every
  breakpoint and stays inside the physical range.
- The shell structure. The level density parameter is suppressed threefold at the doubly magic
  heavy fragment relative to a mid-shell fragment, which is what gives `R_a(A_H)` its structure.
- The shape of the result. The parameterized `r_ν(A_H)` reproduces the published systematic
  behaviour — one half at the symmetric split, a minimum near the heavy magic fragment, one half
  again near the most probable fragmentation, and a near-linear rise above it.
- Sensitivity to the charge polarization. Substituting a tabulated polarization for the average
  values moves `R_T(A_H)` by at most 3.4 %, with a median of 0.27 %, worst at the shell minimum.
- The suite passes on the declared Julia floor and on the current release, 223 assertions on each.
- The back-shifted Fermi gas prescription, against the published text. The three coefficients, the
  shell correction, the deuteron pairing term and all five liquid-drop coefficients reproduce
  Phys. Rev. C **72**, 044311 (2005), Eqs. (7) and (9), and Phys. Rev. C **80**, 054310 (2009),
  Eqs. (12) and (19).

- The shell corrections of the Gilbert-Cameron prescription are read in the order the table
  declares them and looked up at the right nucleon number, `S(Z)` at the proton number and `S(N)`
  at the neutron number. Exchanging the two columns is not absorbed by their sum, and would be
  silent; a test pins the order against the first tabulated row and the closed form.

Not verified:

- **No published number has been reproduced.** The total average temperature ratio quoted in the
  literature is taken over a fission fragment mass yield distribution `Y(A)`, which this package
  does not take as input, so the one directly comparable quantity cannot yet be computed. The
  agreement established so far is of shape, not of value.
- The Gilbert-Cameron prescription has been exercised for magnitude, for its expected departure
  from the back-shifted Fermi gas, and for the table lookup, but its values have not been compared
  against a published tabulation.
- The effect of neglecting the back-shift. The level density parameter is taken from a systematic
  that fits it jointly with a back-shift `E1`, while the extraction rests on the un-shifted
  `E* = a T²` the method is published under. `E1` differs between the two fragments, so it does not
  cancel in the ratio; at order ±1 MeV against fragment excitations of 10-20 MeV the effect is
  expected to be small, but it has not been quantified.
- The continuous integration workflow has never run.
- The path for a fissioning nucleus of odd mass number is covered only by unit tests; no such case
  exists in the data.

## Status

| Component | State |
|---|---|
| Mass excess and charge distribution input | complete |
| Level density parameter, back-shifted Fermi gas | complete |
| Level density parameter, Gilbert-Cameron | complete |
| Fragmentation range and isobaric charge distribution | complete |
| Multiplicity ratio and temperature ratio | complete |
| Segmented parameterization with model selection | complete |
| Pipeline, tabulated output, figures, provenance | complete |
| Per-data-set and systematic-trend parameterizations | complete |
| Averaging over a fragment mass yield distribution | not implemented; needs Y(A) as input |
| Reproduction of published total averages | blocked on the above |

## Licensing

The source code is under the MIT licence in `LICENSE`.

That licence does not extend to the contents of `data/`, none of which originates with this
package. The atomic mass evaluation, the charge distribution systematics, the shell corrections
and every prompt neutron multiplicity measurement are third-party scientific data, redistributed
here for reproducibility under the terms of their own sources. `data/README.md` records what each
file is, where it came from, and how it should be cited; any result derived from a measurement
should cite that measurement.
