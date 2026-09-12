# Nuclear mass excess input.

"""
    MassExcessTable

Mass excesses `D(A, Z)` in MeV, indexed by `(A, Z)`, together with the proton and neutron mass
excesses that the binding energy of any nuclide is measured against.

# Fields

- `D::Dict{Tuple{Int,Int},Float64}`: mass excess in MeV, keyed by `(A, Z)`.
- `Dᵖ::Float64`: proton mass excess in MeV.
- `Dⁿ::Float64`: neutron mass excess in MeV.
- `source::String`: path of the file the table was read from, recorded for provenance.
"""
struct MassExcessTable
    D::Dict{Tuple{Int,Int},Float64}
    Dᵖ::Float64
    Dⁿ::Float64
    source::String
end

"""
    mass_excess(masses, A, Z) -> Union{Float64,Missing}

Mass excess of the nuclide `(A, Z)` in MeV, or `missing` when it is absent from the table.
"""
function mass_excess(masses::MassExcessTable, A::Integer, Z::Integer)
    return get(masses.D, (Int(A), Int(Z)), missing)
end

"""
    read_mass_excess(path) -> MassExcessTable

Read a whitespace-separated mass excess table with the column layout

```
Z  A  symbol  D  σD
```

where `D` and `σD` are in keV, as distributed in the atomic mass evaluation. Values are converted
to MeV on load, since every downstream energy in this package is in MeV.

Throws an `ArgumentError` naming the file when it is absent, unreadable as this layout, or missing
either the proton or the neutron entry, both of which are required to form binding energies.
"""
function read_mass_excess(path::AbstractString)
    isfile(path) || throw(ArgumentError("mass excess file not found: $(path)"))

    table = try
        CSV.read(
            path,
            DataFrame;
            delim = ' ',
            ignorerepeated = true,
            header = ["Z", "A", "symbol", "D", "σD"],
            types = Dict(:Z => Int, :A => Int, :D => Float64, :σD => Float64),
        )
    catch err
        throw(ArgumentError("mass excess file $(path) does not have the layout \
                             `Z A symbol D σD`: $(err)"))
    end

    D = Dict{Tuple{Int,Int},Float64}()
    for row in eachrow(table)
        D[(row.A, row.Z)] = 1e-3 * row.D
    end

    haskey(D, (1, 1)) ||
        throw(ArgumentError("mass excess file $(path) has no proton entry (A = 1, Z = 1)"))
    haskey(D, (1, 0)) ||
        throw(ArgumentError("mass excess file $(path) has no neutron entry (A = 1, Z = 0)"))

    return MassExcessTable(D, D[(1, 1)], D[(1, 0)], String(path))
end
