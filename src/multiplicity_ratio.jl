# Experimental prompt neutron multiplicity data and the multiplicity ratio derived from it.

"""
    Multiplicity

An experimental prompt neutron multiplicity distribution `ν(A)` with its uncertainties.

# Fields

- `A`: fragment mass numbers, ascending.
- `ν`: prompt neutron multiplicity.
- `σν`: uncertainty of `ν`; `missing` where the source quotes none.
- `label`: identifier of the dataset, used in figure legends and output file names.
- `source`: path of the file the data was read from, recorded for provenance.
"""
struct Multiplicity
    A::Vector{Int}
    ν::Vector{Float64}
    σν::Vector{Union{Missing,Float64}}
    label::String
    source::String
end

Base.length(data::Multiplicity) = length(data.A)

"""
    multiplicity(data, A) -> Union{Tuple{Float64,Float64},Missing}

Multiplicity and its uncertainty at mass number `A`, or `missing` if the dataset has no entry
there.

# Examples

```jldoctest
julia> data = Multiplicity([120, 132], [3.10, 0.69], [0.05, 0.03], "example", "");

julia> multiplicity(data, 132)
(0.69, 0.03)

julia> multiplicity(data, 131)
missing
```
"""
function multiplicity(data::Multiplicity, A::Integer)
    index = findfirst(==(A), data.A)
    return index === nothing ? missing : (data.ν[index], data.σν[index])
end

"""
    read_multiplicity(path; label) -> Multiplicity

Read a whitespace-separated `ν(A)` dataset with the column layout

```
A  nu  nu_uncertainty
```

and a single header line, or `A nu` where the source quotes no uncertainty. An absent third
field, and one that is zero, negative or not finite (`NaN`), is stored as `missing`: an
unquoted uncertainty is not a datum, and it is never written as zero, which would denote an exact
value.

The columns are taken **by position**, not by header text: the header line is skipped, so the
upstream retrieval may rename it without touching anything here, and nothing in this reader may
be changed to a lookup by name. Fields beyond the third are not read, and blank lines are passed
over.

No point is ever dropped on the basis of its value or uncertainty: filtering experimental data
requires a documented reason specific to the dataset, which belongs with the data rather than in
this reader.

Throws an `ArgumentError` naming the file, and the line where one is at fault, when a line holds
fewer than two fields, when a mass number, a multiplicity or an uncertainty is not numeric, when a
multiplicity is negative or not finite, when mass numbers repeat, or when the file holds no data.
A field that cannot be read is refused rather than taken as unquoted: one malformed uncertainty
must not pass for a dataset that quotes none.
"""
function read_multiplicity(path::AbstractString; label::AbstractString = "")
    isfile(path) || throw(ArgumentError("multiplicity file not found: $(path)"))

    A = Int[]
    ν = Float64[]
    σν = Union{Missing,Float64}[]
    for (number, line) in enumerate(eachline(path))
        number == 1 && continue
        fields = split(line)
        isempty(fields) && continue
        malformed(what) = throw(
            ArgumentError("multiplicity file $(path), line $(number): $(what); the layout is \
                 `A nu nu_uncertainty`, got $(repr(line))"),
        )
        length(fields) ≥ 2 || malformed("fewer than two fields")
        mass = tryparse(Float64, fields[1])
        (mass === nothing || !isfinite(mass)) &&
            malformed("the mass number $(repr(fields[1])) is not numeric")
        value = tryparse(Float64, fields[2])
        value === nothing && malformed("the multiplicity $(repr(fields[2])) is not numeric")
        (isfinite(value) && value ≥ 0) ||
            malformed("the multiplicity $(repr(fields[2])) is negative or not finite")
        uncertainty = if length(fields) ≥ 3
            quoted = tryparse(Float64, fields[3])
            quoted === nothing &&
                malformed("the uncertainty $(repr(fields[3])) is not numeric")
            # Zero, a negative number and NaN are how a source writes "none quoted".
            (isfinite(quoted) && quoted > 0) ? quoted : missing
        else
            missing
        end
        push!(A, round(Int, mass))
        push!(ν, value)
        push!(σν, uncertainty)
    end

    allunique(A) || throw(ArgumentError("multiplicity file $(path) has repeated mass numbers"))
    isempty(A) && throw(ArgumentError("multiplicity file $(path) contains no usable rows"))

    order = sortperm(A)
    name = isempty(label) ? splitext(basename(path))[1] : String(label)
    return Multiplicity(A[order], ν[order], σν[order], name, String(path))
end

"""
    RatioCurve

A ratio tabulated against the heavy-fragment mass number, with uncertainties.

Used both for the prompt neutron multiplicity ratio `r_ν(A_H)` and for the temperature ratio
`R_T(A_H)`.

# Fields

- `A_H`: heavy-fragment mass numbers, ascending.
- `ratio`: the ratio itself. Named for the quantity rather than for its role, which is what lets
  the tables written from it carry a column named `r_nu` or `R_T` instead of `value`.
- `σ`: uncertainty of the ratio; `missing` where the data it came from quote none, and zero
  only where the value is exact, at a pinned abscissa.
- `label`: identifier inherited from the underlying dataset.
"""
struct RatioCurve
    A_H::Vector{Int}
    ratio::Vector{Float64}
    σ::Vector{Union{Missing,Float64}}
    label::String
end

Base.length(curve::RatioCurve) = length(curve.A_H)
Base.isempty(curve::RatioCurve) = isempty(curve.A_H)

"""
    multiplicity_ratio(data, A₀, A_H_range) -> RatioCurve

The prompt neutron multiplicity ratio

```
r_ν(A_H) = ν_H / (ν_L + ν_H),    ν_H = ν(A_H),  ν_L = ν(A₀ - A_H),
```

for every heavy mass number of `A_H_range` at which the dataset provides both fragments of the
pair. This is the ratio that equals `E*_H / TXE` under the assumption that the multiplicity ratio
of complementary fragments follows their excitation energy ratio.

Uncertainties are propagated from those of `ν`,

```
σ_r² = [ν_L² σ_H² + ν_H² σ_L²] / (ν_L + ν_H)⁴.
```

Where one fragment of the pair quotes no uncertainty its term is dropped and the result is a lower
bound carried as quoted, which is how the published values were obtained from such datasets;
where neither quotes one the result is `missing`.

Pairs whose ratio falls outside the open interval `(0, 1)` are omitted rather than clipped: the
endpoints are the singular points of the temperature ratio relation, so a clipped value would
enter the extraction as a spurious datum.
"""
function multiplicity_ratio(data::Multiplicity, A₀::Integer, A_H_range::UnitRange{Int})
    A_H = Int[]
    ratio = Float64[]
    σ = Union{Missing,Float64}[]

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
        σ_r = if ismissing(σ_H) && ismissing(σ_L)
            missing
        else
            sqrt((ν_L * coalesce(σ_H, 0.0))^2 + (ν_H * coalesce(σ_L, 0.0))^2) / total^2
        end

        push!(A_H, mass)
        push!(ratio, r)
        push!(σ, σ_r)
    end

    return RatioCurve(A_H, ratio, σ, data.label)
end
