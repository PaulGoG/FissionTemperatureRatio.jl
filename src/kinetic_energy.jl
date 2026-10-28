# The mean total kinetic energy against heavy-fragment mass, which weights the charge-resolved
# inversion by the excitation each fragmentation shares out.

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
- `record`: its [`RetrievalRecord`](@ref), or `nothing` where no record beside it lists the file.
"""
struct MeanKineticEnergy
    values::Dict{Int,Float64}
    measured::Vector{Int}
    averaged::Vector{Int}
    interpolated::Vector{Int}
    extrapolated::Vector{Int}
    label::String
    source::String
    record::Union{RetrievalRecord,Nothing}
end

"""
    read_mean_kinetic_energy(path, system, heavy_masses; record) -> MeanKineticEnergy

Read a `⟨TKE⟩(A)` dataset — one header row, then `A TKE [TKE_uncertainty]`, whitespace-separated,
`TKE` in MeV and pre-neutron, as the EXFOR retrieval writes `data/<system>/TKE_vs_A/` — and carry
it onto every heavy mass number of `heavy_masses`; see [`MeanKineticEnergy`](@ref). Columns are
taken by position. `record` is the dataset's retrieval record, by default the one beside the
file, [`retrieval_record`](@ref).

Throws an `ArgumentError` naming the file when a mass number repeats, a value is not a positive
finite energy, or no tabulated mass falls on the heavy wing of `system` directly or through the
complement.
"""
function read_mean_kinetic_energy(
    path::AbstractString,
    system::FissioningSystem,
    heavy_masses::AbstractUnitRange{<:Integer};
    record::Union{RetrievalRecord,Nothing} = retrieval_record(path),
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
        values, measured, averaged, interpolated, extrapolated, label, String(path), record
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
    for A_H in sort!(collect(keys(energies.values)))
        TKE = energies.values[A_H]
        entry = mass_yield(yields, A_H)
        entry === nothing && continue
        numerator += entry[1] * TKE
        total += entry[1]
        push!(masses, A_H)
    end
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
        standard = Measurements.value(standard),
        standard_uncertainty = Measurements.uncertainty(standard),
        offset = mean - Measurements.value(standard),
        masses,
    )
end
