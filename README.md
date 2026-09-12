# FissionTemperatureRatio.jl

Extraction of the temperature ratio `R_T = T_L/T_H` of complementary fully accelerated fission
fragments from experimental prompt neutron multiplicity data, and its parameterization by joined
straight segments.

```
FissionTemperatureRatio/
├── activate.jl                     activate and instantiate the root environment
├── Project.toml, Manifest.toml     the pinned environment
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

## Licensing

The source code is under the MIT licence in `LICENSE`.

That licence does not extend to the contents of `data/`, none of which originates with this
package. The atomic mass evaluation, the charge distribution systematics, the shell corrections
and every prompt neutron multiplicity measurement are third-party scientific data, redistributed
here for reproducibility under the terms of their own sources. `data/README.md` records what each
file is, where it came from, and how it should be cited; any result derived from a measurement
should cite that measurement.
