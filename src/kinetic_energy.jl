# The mean total kinetic energy against heavy-fragment mass, which weights the charge-resolved
# inversion by the excitation each fragmentation shares out, and the retrieval records that say
# where every input dataset came from.

"""
    FLAGGED_QUALIFIERS

Reaction-code qualifiers that mark a dataset as other than a direct measurement of the quantity
at the entrance channel's energy: `DERIV`, derived from other data, and `SPA`, averaged over an
unspecified neutron spectrum. A dataset carrying either is used, but flagged in the log, in the
dataset diagnostics and in the run metadata.
"""
const FLAGGED_QUALIFIERS = ("DERIV", "SPA")

"""
    RetrievalRecord

The entry of one dataset in the `retrieval.toml` that the EXFOR retrieval writes beside the files
it produced, with the run of that retrieval.

# Fields

- `path`: the record file.
- `qualifiers`: the reaction-code qualifiers the retrieval recorded for the dataset, as tags
  (`"DERIV"`, `"MXW"`, …).
- `entry`: the dataset's `[[accepted]]` table, as written.
- `run`: the record's `[run]` table: configuration, parser version and revision, timestamp.
  Empty for a record written before the parser recorded it.
"""
struct RetrievalRecord
    path::String
    qualifiers::Vector{String}
    entry::Dict{String,Any}
    run::Dict{String,Any}
end

# Every `retrieval*.toml` in a directory: one per query that contributed files to it.
function _retrieval_documents(directory::AbstractString)
    isdir(directory) || return Pair{String,Dict{String,Any}}[]
    names = sort!(
        filter(f -> startswith(f, "retrieval") && endswith(f, ".toml"), readdir(directory))
    )
    return [
        joinpath(directory, name) => TOML.parsefile(joinpath(directory, name)) for name in names
    ]
end

# "MXW: Maxwellian-averaged" → "MXW"; the retrieval writes the tag and its meaning.
_qualifier_tag(qualifier::AbstractString) = String(strip(first(split(qualifier, ':'))))

"""
    retrieval_record(path) -> Union{RetrievalRecord,Nothing}

The retrieval record of the data file at `path`: its `[[accepted]]` entry in a `retrieval*.toml`
of the same directory, or `nothing` where no record there lists the file.
"""
function retrieval_record(path::AbstractString)
    file = basename(path)
    for (record, document) in _retrieval_documents(dirname(path))
        for entry in get(document, "accepted", Any[])
            entry isa AbstractDict && get(entry, "file", nothing) == file || continue
            qualifiers = [_qualifier_tag(q) for q in get(entry, "qualifiers", String[])]
            run = get(document, "run", Dict{String,Any}())
            return RetrievalRecord(
                record, qualifiers, Dict{String,Any}(entry), Dict{String,Any}(run)
            )
        end
    end
    return nothing
end

"""
    retrieval_qualifiers(directory) -> Dict{String,Vector{String}}

The reaction-code qualifier tags of every data file the retrieval records in `directory` list,
keyed by file name. A file no record lists is absent.
"""
function retrieval_qualifiers(directory::AbstractString)
    qualifiers = Dict{String,Vector{String}}()
    for (_, document) in _retrieval_documents(directory)
        for entry in get(document, "accepted", Any[])
            entry isa AbstractDict && haskey(entry, "file") || continue
            qualifiers[entry["file"]] = [
                _qualifier_tag(q) for q in get(entry, "qualifiers", String[])
            ]
        end
    end
    return qualifiers
end

"""
    pooling_weight(record) -> Float64

The factor by which a dataset's values are weighted where datasets are pooled: raw points over
rows written, `mass_values_non_integer / rows_written`, for a dataset whose retrieval interpolated
non-integer masses onto the integers (`mass_treatment = "interpolated"`), and one otherwise or
where no record is held. Neighbouring interpolated rows share their bracketing points, so a dataset
written at more masses than it measured carries no more information than its measurements.
"""
pooling_weight(::Nothing) = 1.0
function pooling_weight(record::RetrievalRecord)
    get(record.entry, "mass_treatment", nothing) == "interpolated" || return 1.0
    points = get(record.entry, "mass_values_non_integer", nothing)
    rows = get(record.entry, "rows_written", nothing)
    (points isa Integer && rows isa Integer && 0 < points <= rows) || return 1.0
    return points / rows
end

"""
    pair_sum_scale(record) -> Union{NamedTuple,Nothing}

The scale of a multiplicity dataset as the retrieval states it, for one read by the complement
test: `deviation`, the relative departure of its pair sum `ν(A) + ν(A₀ - A)`, weighted with the
light-fragment yield, from `ν̄`; its `uncertainty`; and `consistent`, whether the two agree within
three standard deviations. `uncertainty` and `consistent` are `missing` where the dataset states
no uncertainty, and the result is `nothing` where the record forms no pair sum.

A uniform scale cancels in `r_ν = ν_H/(ν_L + ν_H)`, so a dataset off the scale of `ν̄` is used
and flagged, not corrected. A scale error that varies with mass does not cancel, and the pair sum
cannot detect it; the structural diagnostics of the dataset are what bear on it.
"""
pair_sum_scale(::Nothing) = nothing
function pair_sum_scale(record::RetrievalRecord)
    deviation = get(record.entry, "pair_sum_deviation", nothing)
    deviation isa Real || return nothing
    uncertainty = get(record.entry, "pair_sum_deviation_uncertainty", missing)
    consistent = get(record.entry, "scale_consistent", missing)
    return (
        deviation = Float64(deviation),
        uncertainty = uncertainty isa Real ? Float64(uncertainty) : missing,
        consistent = consistent isa Bool ? consistent : missing,
    )
end

function _flagged(qualifiers::AbstractVector{<:AbstractString})
    return [q for q in qualifiers if q in FLAGGED_QUALIFIERS]
end

"""
    MeanKineticEnergy

The pre-neutron mean total kinetic energy `⟨TKE⟩(A_H)` in MeV at every heavy mass number of a
fragmentation range, from one measured dataset.

`TKE` is a property of the pair: the two fragments of a split are counted in one event, so
`⟨TKE⟩(A) = ⟨TKE⟩(A₀ - A)` exactly for pre-neutron masses. Where the dataset gives both `A_H` and
its complement, the heavy mass takes their mean, so that a difference between the two wings — a
backing loss on one side of a double-energy measurement — enters once and alike; where it gives
one, that value serves the pair. A heavy mass between two measured ones takes the linear
interpolation, and one beyond the measured span takes the value at the nearest measured mass.
Each is listed, so that a run can say which of its weights rest on data at their own mass number.

# Fields

- `values`: `⟨TKE⟩` in MeV, keyed by heavy mass number, at every mass of the range.
- `measured`: heavy mass numbers of the range the dataset gives, directly or through the
  complement.
- `averaged`: those of them the dataset gives on both wings, whose value is the mean of the two.
- `interpolated`, `extrapolated`: heavy mass numbers of the range filled from neighbours, and
  from the nearest measured mass beyond the measured span.
- `label`: the dataset, named as the multiplicity datasets are.
- `source`: the file read.
"""
struct MeanKineticEnergy
    values::Dict{Int,Float64}
    measured::Vector{Int}
    averaged::Vector{Int}
    interpolated::Vector{Int}
    extrapolated::Vector{Int}
    label::String
    source::String
end

"""
    read_mean_kinetic_energy(path, system, heavy_masses) -> MeanKineticEnergy

Read a `⟨TKE⟩(A)` dataset — one header row, then `A TKE [TKE_uncertainty]`, whitespace-separated,
`TKE` in MeV and pre-neutron, as the EXFOR retrieval writes `data/<system>/TKE_vs_A/` — and carry
it onto every heavy mass number of `heavy_masses`; see [`MeanKineticEnergy`](@ref). Columns are
taken by position.

Throws an `ArgumentError` naming the file when a mass number repeats, a value is not a positive
finite energy, or no tabulated mass falls on the heavy wing of `system` directly or through the
complement.
"""
function read_mean_kinetic_energy(
    path::AbstractString, system::FissioningSystem, heavy_masses::AbstractUnitRange{<:Integer}
)
    table = read_delimited_table(path, MEAN_KINETIC_ENERGY_SPEC)
    masses = integer_column(table, :A)
    energies = column(table, :TKE, Float64)
    allunique(masses) || throw(ArgumentError("$(path) repeats a mass number"))
    for (A, TKE) in zip(masses, energies)
        isfinite(TKE) && TKE > 0 || throw(
            ArgumentError(
                "$(path) has a non-positive or non-finite ⟨TKE⟩ = $(TKE) at A = $(A)"
            ),
        )
    end

    A₀ = system.compound.A
    heavy = Dict{Int,Float64}()
    complement = Dict{Int,Float64}()
    for (A, TKE) in zip(masses, energies)
        2 * A >= A₀ ? (heavy[A] = TKE) : (complement[A₀ - A] = TKE)
    end
    # Both wings of one split: their mean. One wing: it serves the pair.
    both = sort!(collect(intersect(keys(heavy), keys(complement))))
    tabulated = merge(complement, heavy)
    for A_H in both
        tabulated[A_H] = (heavy[A_H] + complement[A_H]) / 2
    end
    isempty(tabulated) && throw(ArgumentError("$(path) tabulates no mass number"))
    known = sort!(collect(keys(tabulated)))

    values = Dict{Int,Float64}()
    measured = Int[]
    interpolated = Int[]
    extrapolated = Int[]
    for A_H in heavy_masses
        if haskey(tabulated, A_H)
            values[A_H] = tabulated[A_H]
            push!(measured, A_H)
        elseif A_H < first(known)
            values[A_H] = tabulated[first(known)]
            push!(extrapolated, A_H)
        elseif A_H > last(known)
            values[A_H] = tabulated[last(known)]
            push!(extrapolated, A_H)
        else
            upper = known[searchsortedfirst(known, A_H)]
            lower = known[searchsortedlast(known, A_H)]
            fraction = (A_H - lower) / (upper - lower)
            values[A_H] = (1 - fraction) * tabulated[lower] + fraction * tabulated[upper]
            push!(interpolated, A_H)
        end
    end
    isempty(measured) && throw(
        ArgumentError(
            "$(path) gives ⟨TKE⟩ at no heavy mass number of \
             $(first(heavy_masses)):$(last(heavy_masses)), directly or through the complement",
        ),
    )
    label = _dataset_label(basename(path))
    averaged = [A_H for A_H in measured if A_H in both]
    return MeanKineticEnergy(
        values, measured, averaged, interpolated, extrapolated, label, String(path)
    )
end

"""
    mean_kinetic_energy_offset(energies, yields, system) -> NamedTuple

The yield-weighted mean of a `⟨TKE⟩(A_H)` dataset over the heavy masses it was carried to,
`Σ Y(A_H) ⟨TKE⟩(A_H) / Σ Y(A_H)` with the mass yields `yields`, and its offset from the energy
standard FissionFragmentsDomain records for `system` (Gönnenwein's recommendations as tabulated by
Bertsch et al., J. Phys. G 42, 077001 (2015), doi:10.1088/0954-3899/42/7/077001). Double-energy
measurements differ in their pulse-height-defect calibration by several MeV, so the offset says
how far the input sits from the scale the literature normalises to, provided `yields` spans the
heavy peak: a distribution measured over part of it weights the mean towards its own masses. All
energies in MeV: `(mean, standard, standard_uncertainty, offset, masses)`, the standard, its
uncertainty and the offset `nothing` where no standard exists, `masses` the heavy mass numbers
`yields` gives, over which the mean was taken. Throws where `yields` covers none of the masses.
"""
function mean_kinetic_energy_offset(
    energies::MeanKineticEnergy, yields::MassYield, system::FissioningSystem
)
    numerator = 0.0
    total = 0.0
    masses = Int[]
    for (A_H, TKE) in energies.values
        entry = mass_yield(yields, A_H)
        entry === nothing && continue
        numerator += entry[1] * TKE
        total += entry[1]
        push!(masses, A_H)
    end
    sort!(masses)
    total > 0 || throw(
        ArgumentError("the yield distribution $(repr(yields.label)) covers no heavy mass of \
             $(repr(energies.label))"),
    )
    mean = numerator / total
    standard = recommended_mean_total_kinetic_energy(system)
    standard === nothing && return (;
        mean, standard = nothing, standard_uncertainty = nothing, offset = nothing, masses
    )
    return (;
        mean,
        standard = value(standard),
        standard_uncertainty = uncertainty(standard),
        offset = mean - value(standard),
        masses,
    )
end
