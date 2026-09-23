# Fission fragment mass yields, and the yield-weighted total average of a ratio curve.

"""
    MassYield

An experimental fission fragment mass yield distribution `Y(A)` with its uncertainties.

The yields are pre-neutron: the temperature ratio is a function of the primary heavy-fragment mass
number, so averaging it over a post-neutron distribution would weight each ratio by the yield of a
different fragmentation.

No normalization is imposed. A total average divides by the sum of the weights it used, so a
distribution summing to unity over one peak, to two over the full mass range, or to a hundred per
cent gives the same answer.

# Fields

- `A`: fragment mass numbers, ascending.
- `Y`: mass yield.
- `σY`: uncertainty of `Y`; `missing` where the source quotes none.
- `label`: identifier of the distribution, used in output tables and figure legends.
- `source`: path of the file the data was read from, recorded for provenance.
"""
struct MassYield
    A::Vector{Int}
    Y::Vector{Float64}
    σY::Vector{Union{Missing,Float64}}
    label::String
    source::String
end

Base.length(data::MassYield) = length(data.A)
Base.isempty(data::MassYield) = isempty(data.A)

"""
    mass_yield(data, A) -> Union{Tuple{Float64,Float64},Missing}

Yield and its uncertainty at mass number `A`, or `missing` if the distribution has no entry there.

# Examples

```jldoctest
julia> data = MassYield([132, 140], [0.061, 0.048], [0.002, missing], "example", "");

julia> mass_yield(data, 132)
(0.061, 0.002)

julia> mass_yield(data, 133)
missing
```
"""
function mass_yield(data::MassYield, A::Integer)
    index = findfirst(==(A), data.A)
    return index === nothing ? missing : (data.Y[index], data.σY[index])
end

"""
    read_mass_yield(path; label) -> MassYield

Read a whitespace-separated `Y(A)` dataset with the column layout

```
A  Y  Y_uncertainty
```

and a single header line, or `A Y` where the source quotes no uncertainty. An absent,
non-numeric or non-positive third column is stored as `missing`, never as zero.

The columns are taken **by position**, not by header text: the header line is skipped, so the
upstream retrieval may rename it without touching anything here, and nothing in this reader may
be changed to a lookup by name.

As for multiplicity data, no point is dropped on the basis of its value. A negative yield is an
error rather than a datum, and is rejected.

Throws an `ArgumentError` naming the file when it cannot be read with this layout, when mass
numbers repeat, or when a yield is negative.
"""
function read_mass_yield(path::AbstractString; label::AbstractString = "")
    isfile(path) || throw(ArgumentError("yield file not found: $(path)"))

    table = try
        CSV.read(
            path,
            DataFrame;
            delim = ' ',
            ignorerepeated = true,
            header = ["A", "Y", "σY"],
            skipto = 2,
            silencewarnings = true,
        )
    catch exception
        throw(ArgumentError("yield file $(path) does not have the layout \
                             `A Y Y_uncertainty`: $(exception)"))
    end

    A = Int[]
    Y = Float64[]
    σY = Union{Missing,Float64}[]
    for row in eachrow(table)
        (ismissing(row.A) || ismissing(row.Y)) && continue
        mass = round(Int, row.A)
        value = Float64(row.Y)
        value ≥ 0 ||
            throw(ArgumentError("yield file $(path) has a negative yield at A = $(mass)"))
        push!(A, mass)
        push!(Y, value)
        push!(σY, _quoted_uncertainty(table, row, :σY))
    end

    allunique(A) || throw(ArgumentError("yield file $(path) has repeated mass numbers"))
    isempty(A) && throw(ArgumentError("yield file $(path) contains no usable rows"))

    order = sortperm(A)
    name = isempty(label) ? splitext(basename(path))[1] : String(label)
    return MassYield(A[order], Y[order], σY[order], name, String(path))
end

"""
    read_mass_yield_directory(directory) -> Vector{MassYield}

Read every `.dat` file in `directory`, sorted by name so that output rows are ordered
reproducibly. Other files are ignored, so a run record can sit beside the data it describes.

The label of each distribution is its file name stripped of the leading archive identifier and the
extension, with underscores replaced by spaces.
"""
function read_mass_yield_directory(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("yield directory not found: $(directory)"))
    files = _data_files(directory)
    isempty(files) && throw(ArgumentError("yield directory holds no data files: $(directory)"))
    return [
        read_mass_yield(joinpath(directory, file); label = _dataset_label(file)) for
        file in files
    ]
end

"""
    total_average(curve, yields) -> Tuple{Float64,Float64}

The total average of a ratio curve over a fragment mass yield distribution,

```
⟨R_T⟩ = Σ Y(A_H) R_T(A_H) / Σ Y(A_H),
```

with its uncertainty, taken over the heavy mass numbers the two have in common.

This is the quantity the literature tabulates, and the one a prompt emission code takes when it
uses a single temperature ratio for all fragmentations rather than a function of mass number. It
differs from the mean over the fragment mass range, which weights every mass number equally and
so is dominated by the far-asymmetric tail where the yield is negligible.

The normalization of `Y` cancels, so a distribution summing to unity, two, or a hundred gives the
same result.

Both inputs carry uncertainties and both propagate:

```
σ² = Σ [ (Y_i/ΣY)² σ_{R,i}² + ((R_i - ⟨R_T⟩)/ΣY)² σ_{Y,i}² ]
```

The second term is what makes the average of a multiplicity dataset quoting no uncertainties
still carry one, from the yield distribution alone. An unquoted uncertainty, `missing` in either
input, contributes nothing; where neither input quotes any the uncertainty is zero. The yield
term is a difference from the mean, so a distribution contributes nothing to the uncertainty at
mass numbers where the ratio sits at its own average.

Correlations between mass numbers are neglected in both inputs. For the ratio that is the
independent-points approximation of the published tables, and it is exact only for a curve whose
points are measured independently; the tabulated values of a segmented curve are functions of a
few fitted coefficients, and the method for a [`SegmentedCurve`](@ref) propagates their
covariance instead.

Throws an `ArgumentError` when the curve and the distribution share no mass number, or when the
yields summed over the shared mass numbers are not positive.
"""
function total_average(curve::RatioCurve, yields::MassYield)
    weight = Float64[]
    ratio = Float64[]
    σ_R = Float64[]
    σ_Y = Float64[]

    for (index, mass) in enumerate(curve.A_H)
        entry = mass_yield(yields, mass)
        ismissing(entry) && continue
        push!(weight, entry[1])
        push!(σ_Y, coalesce(entry[2], 0.0))
        push!(ratio, curve.ratio[index])
        push!(σ_R, coalesce(curve.σ[index], 0.0))
    end

    isempty(weight) &&
        throw(ArgumentError("the ratio curve \"$(curve.label)\" and the yield distribution \
             \"$(yields.label)\" share no mass number"))
    total = sum(weight)
    total > 0 || throw(
        ArgumentError(
            "the yield distribution \"$(yields.label)\" sums to $(total) over the mass \
             numbers it shares with \"$(curve.label)\""
        ),
    )

    mean = dot(weight, ratio) / total
    variance = 0.0
    for index in eachindex(weight)
        variance += (weight[index] / total)^2 * σ_R[index]^2
        variance += ((ratio[index] - mean) / total)^2 * σ_Y[index]^2
    end
    return (mean, sqrt(variance))
end
