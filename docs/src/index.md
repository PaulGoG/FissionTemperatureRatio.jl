# FissionTemperatureRatio.jl

FissionTemperatureRatio.jl extracts the temperature ratio ``R_T = T_L/T_H`` of complementary,
fully accelerated fission fragments from measured prompt neutron multiplicities ``\nu(A)`` and
describes it by joined straight segments.

Prompt emission model codes that share the total excitation energy at full acceleration take
``R_T`` as input. Obtaining it by fitting ``\nu(A)`` with such a code makes the result dependent
on the emission model of that code; the method implemented here uses the same experimental
``\nu(A)`` and no emission calculation, so its result is independent of any particular prompt
emission treatment.

The fragmentation domain, the level density parameters and the relation between ``R_T`` and the
excitation-energy partition are those of
[FissionFragmentsDomain.jl](https://PaulGoG.github.io/FissionFragmentsDomain.jl/stable/), shared
with the prompt emission codes that read ``R_T(A_H)``, so that the curve is extracted on the
domain it is applied on. The [method](method.md) page sets out how.

## Installation

The package is not registered. Clone it; every environment activates and instantiates itself:

```
julia -i activate.jl
```

## Running the pipeline

```julia
using FissionTemperatureRatio

configuration = load_configuration("config/U233_nth.toml"; data_directory = "data")
result = run_pipeline(configuration)
write_results(result, joinpath("data", "sims", "U233_nth", run_identifier(configuration)))
```

or, from a shell,

```
julia scripts/run.jl config/U233_nth.toml
```

[`run_pipeline`](@ref) returns an [`ExtractionResult`](@ref) and writes nothing.
[`write_results`](@ref) writes the tables and the manifest into the directory it is given, and
refuses a directory that already holds files.

The script writes to `data/sims/<system>/<run identifier>/`. If a run with the same identifier
exists, it is first moved aside as `<run identifier>#1`, `#2`, …. Beside the tables the script
writes the provenance record `metadata.toml`, a copy of the configuration as `configuration.toml`
and a copy of the resolved manifest of its environment as `Manifest.toml`. The figures go to
`plots/<system>/<run identifier>/`. The run identifier covers every configuration key that
changes the result, so the directory name identifies the run. Four files inside it repeat the
identifier in their names: `manifest_<run identifier>.toml`,
`segmented_curves_<run identifier>.csv`, `total_average_R_T_<run identifier>.csv` and
`leave_one_out_<run identifier>.csv`.

A consuming code stages the whole directory and selects the manifest by its prefix. The manifest
records the domain the curves were extracted on, and a consuming code refuses a curve extracted
on a different one. [Naming](naming.md) describes the file layout, and
[Configuration](configuration.md) lists every configuration key with its type, bounds, default
and run-identifier token.
