# Prototype

This branch is the state of the work before it was rewritten as a package. It is kept for the
record and is not developed: nothing here is maintained, and no change should be committed to it.
The rewrite lives on `main`.

```
functionBodies.jl   level density parameter, fragmentation domain, ratios, the segment fit
preamble.jl         hardwired parameters: nucleus, mass range, charges per mass, fit settings
main.jl             the run, top to bottom
```

It works, and it produced results. What it is not: there is no project environment and no
manifest, so the dependency versions it ran against are unrecorded; there are no tests; every
parameter is a `const` in the source, so a different fissioning nucleus means editing the code;
paths are relative to a `cd` at load; and results are written over their predecessors with no run
identifier or commit recorded. Several of its numerical choices also depart from the published
method — the averaging order of the level density parameter ratio, the segment-selection
objective, the treatment of uncertainties through the fit, and an undocumented filter on the input
data. `CHANGELOG.md` on `main` lists them under what the rewrite corrected.

## Input data

The input data is no longer tracked here. It remains reachable in this branch's history, and in
the commits `main` shares with it, so nothing is lost; it is simply not part of what the branch
ships. `data/README.md` on `main` records what each input is, where it comes from and under what
terms, and the experimental measurements are now retrieved with
[ExforFissionData.jl](https://github.com/PaulGoG/ExforFissionData.jl).

To run this code you would need to restore `inputData/` from history:

```
git checkout 1c5d3b5 -- inputData
```

## Method

A. Tudora, P. Gogita, *Eur. Phys. J. A* **60**, 190 (2024),
[doi:10.1140/epja/s10050-024-01375-7](https://doi.org/10.1140/epja/s10050-024-01375-7).
