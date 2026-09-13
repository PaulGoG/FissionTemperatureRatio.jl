# Nuclear mass excess input.

"""
    MassExcessTable

Mass excesses `Δ(A, Z)` in MeV, indexed by `(A, Z)`, together with the proton and neutron mass
excesses that the binding energy of any nuclide is measured against.

# Fields

- `Δ::Dict{Tuple{Int,Int},Float64}`: mass excess in MeV, keyed by `(A, Z)`.
- `Δᵖ::Float64`: proton mass excess in MeV.
- `Δⁿ::Float64`: neutron mass excess in MeV.
- `source::String`: path of the file the table was read from, recorded for provenance.
"""
struct MassExcessTable
    Δ::Dict{Tuple{Int,Int},Float64}
    Δᵖ::Float64
    Δⁿ::Float64
    source::String
end

"""
    mass_excess(masses, A, Z) -> Union{Float64,Missing}

Mass excess of the nuclide `(A, Z)` in MeV, or `missing` when it is absent from the table.
"""
function mass_excess(masses::MassExcessTable, A::Integer, Z::Integer)
    return get(masses.Δ, (Int(A), Int(Z)), missing)
end

"""
    read_mass_excess_table(path) -> MassExcessTable

Read a whitespace-separated mass excess table with the column layout

```
Z  A  symbol  mass_excess  mass_excess_uncertainty
```

where the last two columns are in keV, as distributed in the atomic mass evaluation. Values are
converted to MeV on load, since every downstream energy in this package is in MeV.

The columns are taken **by position**, not by header text, and this file has no header line at
all: it is the evaluation as distributed. Nothing here may be changed to a lookup by name, which
would couple the reader to a spelling the evaluation never promised.

Throws an `ArgumentError` naming the file when it is absent, unreadable as this layout, or missing
either the proton or the neutron entry, both of which are required to form binding energies.
"""
function read_mass_excess_table(path::AbstractString)
    isfile(path) || throw(ArgumentError("mass excess file not found: $(path)"))

    table = try
        CSV.read(
            path,
            DataFrame;
            delim = ' ',
            ignorerepeated = true,
            header = ["Z", "A", "symbol", "Δ", "σΔ"],
            types = Dict(:Z => Int, :A => Int, :Δ => Float64, :σΔ => Float64),
        )
    catch err
        throw(ArgumentError("mass excess file $(path) does not have the layout \
                             `Z A symbol mass_excess mass_excess_uncertainty`: $(err)"))
    end

    Δ = Dict{Tuple{Int,Int},Float64}()
    for row in eachrow(table)
        Δ[(row.A, row.Z)] = 1e-3 * row.Δ
    end

    haskey(Δ, (1, 1)) ||
        throw(ArgumentError("mass excess file $(path) has no proton entry (A = 1, Z = 1)"))
    haskey(Δ, (1, 0)) ||
        throw(ArgumentError("mass excess file $(path) has no neutron entry (A = 1, Z = 0)"))

    return MassExcessTable(Δ, Δ[(1, 1)], Δ[(1, 0)], String(path))
end

"""
    ShellCorrectionTable

Shell corrections `S_N` and `S_Z` tabulated against nucleon number, as used by the
Gilbert-Cameron level density systematic.

The values are those of Gilbert and Cameron, Can. J. Phys. **43**, 1446 (1965), as distributed in
the IAEA Reference Input Parameter Library, segment on level densities.

# Fields

- `S_N`, `S_Z`: shell correction in MeV, keyed by neutron and proton number respectively.
- `source::String`: path of the file the table was read from, recorded for provenance.
"""
struct ShellCorrectionTable
    S_N::Dict{Int,Float64}
    S_Z::Dict{Int,Float64}
    source::String
end

"""
    read_shell_correction_table(path) -> ShellCorrectionTable

Read a whitespace-separated shell correction table with the column layout

```
n  S_N  S_Z
```

and a single header line, where `n` is read once as a neutron number and once as a proton number.

The columns are taken **by position**, not by header text: the header line is skipped, so
renaming it cannot affect what is read, and nothing here may be changed to a lookup by name.
Exchanging the two correction columns is *not* absorbed by their sum, so the order above is part
of the contract.

Throws an `ArgumentError` naming the file when it is absent or cannot be read with this layout.
"""
function read_shell_correction_table(path::AbstractString)
    isfile(path) || throw(ArgumentError("shell correction file not found: $(path)"))

    table = try
        CSV.read(
            path,
            DataFrame;
            delim = ' ',
            ignorerepeated = true,
            header = ["n", "S_N", "S_Z"],
            skipto = 2,
            types = Dict(:n => Int, :S_N => Float64, :S_Z => Float64),
        )
    catch err
        throw(ArgumentError("shell correction file $(path) does not have the layout \
                             `n S_N S_Z`: $(err)"))
    end

    S_N = Dict{Int,Float64}()
    S_Z = Dict{Int,Float64}()
    for row in eachrow(table)
        S_N[row.n] = row.S_N
        S_Z[row.n] = row.S_Z
    end
    isempty(S_N) &&
        throw(ArgumentError("shell correction file $(path) contains no usable rows"))

    return ShellCorrectionTable(S_N, S_Z, String(path))
end
