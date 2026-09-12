# Experimental prompt neutron multiplicity data and the multiplicity ratio derived from it.

"""
    MultiplicityData

An experimental prompt neutron multiplicity distribution `ν(A)` with its uncertainties.

# Fields

- `A`: fragment mass numbers, ascending.
- `ν`: prompt neutron multiplicity.
- `σν`: uncertainty of `ν`; zero where the source quotes none.
- `label`: identifier of the data set, used in figure legends and output file names.
- `source`: path of the file the data was read from, recorded for provenance.
"""
struct MultiplicityData
    A::Vector{Int}
    ν::Vector{Float64}
    σν::Vector{Float64}
    label::String
    source::String
end

Base.length(data::MultiplicityData) = length(data.A)

"""
    multiplicity(data, A) -> Union{Tuple{Float64,Float64},Missing}

Multiplicity and its uncertainty at mass number `A`, or `missing` if the set has no entry there.
"""
function multiplicity(data::MultiplicityData, A::Integer)
    index = findfirst(==(A), data.A)
    return index === nothing ? missing : (data.ν[index], data.σν[index])
end

"""
    read_multiplicity(path; label) -> MultiplicityData

Read a whitespace-separated `ν(A)` data set with the column layout

```
A  ν  σν
```

and a single header line. A missing or non-numeric third column is taken as an absent
uncertainty and stored as zero, which excludes the point from the weighting of a fit without
discarding it from the plot.

No point is ever dropped on the basis of its value or uncertainty: filtering experimental data
requires a documented reason specific to the data set, which belongs with the data rather than in
this reader.

Throws an `ArgumentError` naming the file when it cannot be read with this layout, when mass
numbers repeat, or when a multiplicity is negative.
"""
function read_multiplicity(path::AbstractString; label::AbstractString = "")
    isfile(path) || throw(ArgumentError("multiplicity file not found: $(path)"))

    table = try
        CSV.read(
            path,
            DataFrame;
            delim = ' ',
            ignorerepeated = true,
            header = ["A", "ν", "σν"],
            skipto = 2,
            silencewarnings = true,
        )
    catch err
        throw(
            ArgumentError(
                "multiplicity file $(path) does not have the layout `A ν σν`: $(err)"
            ),
        )
    end

    A = Int[]
    ν = Float64[]
    σν = Float64[]
    for row in eachrow(table)
        (ismissing(row.A) || ismissing(row.ν)) && continue
        mass = round(Int, row.A)
        value = Float64(row.ν)
        value ≥ 0 || throw(
            ArgumentError(
                "multiplicity file $(path) has a negative multiplicity at A = $(mass)"
            ),
        )
        uncertainty = if hasproperty(table, :σν) && !ismissing(row.σν) && row.σν isa Real
            Float64(row.σν)
        else
            0.0
        end
        push!(A, mass)
        push!(ν, value)
        push!(σν, abs(uncertainty))
    end

    allunique(A) || throw(ArgumentError("multiplicity file $(path) has repeated mass numbers"))
    isempty(A) && throw(ArgumentError("multiplicity file $(path) contains no usable rows"))

    order = sortperm(A)
    name = isempty(label) ? splitext(basename(path))[1] : String(label)
    return MultiplicityData(A[order], ν[order], σν[order], name, String(path))
end

"""
    RatioCurve

A ratio tabulated against the heavy-fragment mass number, with uncertainties.

Used both for the prompt neutron multiplicity ratio `r_ν(A_H)` and for the temperature ratio
`R_T(A_H)`.

# Fields

- `A_H`: heavy-fragment mass numbers, ascending.
- `value`: the ratio.
- `σ`: uncertainty of the ratio.
- `label`: identifier inherited from the underlying data set.
"""
struct RatioCurve
    A_H::Vector{Int}
    value::Vector{Float64}
    σ::Vector{Float64}
    label::String
end

"""
    TREND_LABEL

Label carried by the parameterization that follows the systematic behaviour of the multiplicity
ratio rather than any single experimental data set.
"""
const TREND_LABEL = "systematic trend"

Base.length(curve::RatioCurve) = length(curve.A_H)
Base.isempty(curve::RatioCurve) = isempty(curve.A_H)

"""
    multiplicity_ratio(data, A₀, A_H_range) -> RatioCurve

The prompt neutron multiplicity ratio

```
r_ν(A_H) = ν_H / (ν_L + ν_H),    ν_H = ν(A_H),  ν_L = ν(A₀ - A_H),
```

for every heavy mass number of `A_H_range` at which the data set provides both fragments of the
pair. This is the ratio that equals `E*_H / TXE` under the assumption that the multiplicity ratio
of complementary fragments follows their excitation energy ratio.

Uncertainties are propagated from those of `ν`,

```
σ_r² = [ν_L² σ_H² + ν_H² σ_L²] / (ν_L + ν_H)⁴.
```

Pairs whose ratio falls outside the open interval `(0, 1)` are omitted rather than clipped: the
endpoints are the singular points of the temperature ratio relation, so a clipped value would
enter the extraction as a spurious datum.
"""
function multiplicity_ratio(data::MultiplicityData, A₀::Integer, A_H_range::UnitRange{Int})
    A_H = Int[]
    value = Float64[]
    σ = Float64[]

    for mass in A_H_range
        heavy = multiplicity(data, mass)
        light = multiplicity(data, A₀ - mass)
        (ismissing(heavy) || ismissing(light)) && continue
        ν_H, σ_H = heavy
        ν_L, σ_L = light
        total = ν_L + ν_H
        total > 0 || continue

        r = ν_H / total
        (r > 0 && r < 1) || continue
        σ_r = sqrt((ν_L * σ_H)^2 + (ν_H * σ_L)^2) / total^2

        push!(A_H, mass)
        push!(value, r)
        push!(σ, σ_r)
    end

    return RatioCurve(A_H, value, σ, data.label)
end
