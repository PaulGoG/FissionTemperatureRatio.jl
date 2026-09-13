# FissionTemperatureRatio.jl

Extraction of the temperature ratio ``R_T = T_L/T_H`` of complementary fully accelerated fission
fragments from experimental prompt neutron multiplicity data, and its parameterization by joined
straight segments.

Prompt emission model codes that share the total excitation energy at full acceleration take
``R_T`` as input. Obtaining it by fitting ``\nu(A)`` with such a code makes the result dependent
on the emission model of that code; the method implemented here uses the same experimental
``\nu(A)`` and no emission calculation, so its result is independent of any particular prompt
emission treatment.

## Installation

The package is not registered. Clone it and instantiate the environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Running the pipeline

```julia
using FissionTemperatureRatio

configuration = load_configuration("config/U233_nf.toml")
result = run_pipeline(configuration)
```

or, from a shell,

```
julia --project scripts/run.jl config/U233_nf.toml
```

Results, figures and run metadata are written under `results/` and `plots/`, in a subdirectory
named by the configuration.
