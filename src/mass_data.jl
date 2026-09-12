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

"""
    ShellCorrectionTable

Shell corrections `S(N)` and `S(Z)` tabulated against nucleon number, as used by the
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
    read_shell_corrections(path) -> ShellCorrectionTable

Read a whitespace-separated shell correction table with the column layout

```
n  S(N)  S(Z)
```

and a single header line, where `n` is read once as a neutron number and once as a proton number.

Throws an `ArgumentError` naming the file when it is absent or cannot be read with this layout.
"""
function read_shell_corrections(path::AbstractString)
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
                             `n S(N) S(Z)`: $(err)"))
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
