# FissionTemperatureRatio.jl

Extraction of the temperature ratio ``R_T = T_L/T_H`` of complementary fully accelerated fission
fragments from experimental prompt neutron multiplicity data, described by joined straight
segments.

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

configuration = load_configuration("config/U233_nth.toml")
result = run_pipeline(configuration)
write_results(result, joinpath("data", "sims", "U233_nth", run_identifier(configuration)))
```

or, from a shell,

```
julia scripts/run.jl config/U233_nth.toml
```

[`run_pipeline`](@ref) returns an [`ExtractionResult`](@ref) and writes nothing.
[`write_results`](@ref) writes the tables and the manifest into the directory it is given and
refuses one that already holds files. The script names that directory
`data/sims/<system>/<run identifier>/`, moves an existing run of the same identifier aside as
`<run identifier>#1`, `#2`, …, and adds the provenance record `metadata.toml` beside copies of the
configuration, as `configuration.toml`, and of the resolved manifest of its environment, as
`Manifest.toml`. The figures go to `plots/<system>/<run identifier>/`. The run identifier covers
every configuration key that changes the result, so the directory name alone identifies the run;
inside it, only `manifest_<run identifier>.toml` carries the identifier again, and a consuming code
stages the whole directory and selects the manifest by its prefix. The file layout is set out
under [Naming](naming.md).

The vocabulary the package uses for identifiers, configuration keys, files and column headers is
set out under [Naming](naming.md); it is the same vocabulary the retrieval that supplies the
experimental input writes, and the one a consuming code reads.
