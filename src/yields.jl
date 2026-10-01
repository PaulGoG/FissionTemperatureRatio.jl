# Fission fragment mass yields, and the yield-weighted total average of a ratio curve.
#
# The yield distribution itself, `MassYield`, and its reader are FissionFragmentsDomain's. The
# yields are pre-neutron: the temperature ratio is a function of the primary heavy-fragment mass
# number, so averaging it over a post-neutron distribution would weight each ratio by the yield of
# a different fragmentation. No normalization is imposed; a total average divides by the sum of
# the weights it used.

"""
    read_mass_yield_directory(directory) -> Vector{MassYield}

Read every `.dat` file in `directory`, sorted by name so that output rows are ordered
reproducibly. Other files are ignored, so a run record can sit beside the data it describes.

The label of each distribution is its file name stripped of the leading archive identifier and the
extension, with underscores replaced by spaces. Two files of one author and year carry their
archive identifier in parentheses instead, as the multiplicity datasets do: total averages are
keyed by label, and a shared one would let the second distribution replace the first.
"""
function read_mass_yield_directory(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("yield directory not found: $(directory)"))
    files = _data_files(directory)
    isempty(files) && throw(ArgumentError("yield directory holds no data files: $(directory)"))
    return [
        read_mass_yield(joinpath(directory, file); label = label) for
        (file, label) in zip(files, _unique_labels(files))
    ]
end

"""
    symmetrized_mass_yield(yields, compound_mass) -> MassYield

`yields` with the pre-neutron identity `Y(A) = Y(A₀ - A)` imposed, `A₀` being `compound_mass`: the
two fragments of a split are counted in one event, so their masses have one yield. Where both
complements are measured each takes their mean, with the uncertainty of the mean of two
independent values, `√(σ_A² + σ_{A₀-A}²)/2`; an unquoted one contributes nothing to it, and where
neither is quoted the result quotes none. A mass whose complement is not measured stands for it:
the complement is added with the same yield and uncertainty, so a distribution measured on the
light wing alone gives the heavy-fragment yields a total average is taken over. The result is
ascending in `A`. Where one wing alone was measured its split is counted on both, which an average
does not see, since it divides by the weights it used.

FissionFragmentsDomain's `symmetrized_yield` for a joint `Y(A, TKE)` takes the mean of two measured
complements in the same way, but keeps a cell whose complement is not measured as it is.
"""
function symmetrized_mass_yield(yields::MassYield, compound_mass::Integer)
    A₀ = Int(compound_mass)
    index = Dict(A => i for (i, A) in enumerate(yields.A))
    masses = sort!(union(yields.A, [A₀ - A for A in yields.A if A < A₀]))
    Y = Vector{Float64}(undef, length(masses))
    σY = Vector{Union{Missing,Float64}}(undef, length(masses))
    for (k, A) in enumerate(masses)
        i = get(index, A, nothing)
        j = get(index, A₀ - A, nothing)
        if i === nothing || j === nothing || i == j
            measured = something(i, j)
            Y[k] = yields.Y[measured]
            σY[k] = yields.σY[measured]
        else
            Y[k] = (yields.Y[i] + yields.Y[j]) / 2
            σ = (yields.σY[i], yields.σY[j])
            σY[k] = all(ismissing, σ) ? missing : sqrt(sum(abs2, skipmissing(σ))) / 2
        end
    end
    return MassYield(masses, Y, σY, yields.label, yields.source)
end

"""
    mass_yield_coverage(yields, heavy_masses) -> Float64

The fraction of the heavy mass numbers `heavy_masses` at which `yields` gives a value, on the same
footing as the coverage of a multiplicity dataset. A run takes no total average over a
distribution below `min_dataset_coverage`: averaged over a few heavy masses, `⟨R_T⟩` describes
those masses and not the fission yield. Throws an `ArgumentError` for an empty range.
"""
function mass_yield_coverage(yields::MassYield, heavy_masses::AbstractUnitRange{<:Integer})
    isempty(heavy_masses) && throw(ArgumentError("the heavy mass range is empty"))
    return count(A -> mass_yield(yields, A) !== nothing, heavy_masses) / length(heavy_masses)
end

"""
    yield_fraction(curve, yields, heavy_masses) -> Float64

The fraction of the yield of `yields` over the heavy mass numbers `heavy_masses` that falls at the
mass numbers of `curve`: how much of the distribution a total average over that curve takes in.
One for a curve spanning every mass the distribution covers, less for a dataset whose complete
pairs stop short of them. Throws an `ArgumentError` when the distribution carries no positive
yield in the range.
"""
function yield_fraction(
    curve::RatioCurve, yields::MassYield, heavy_masses::AbstractUnitRange{<:Integer}
)
    on_curve = Set(curve.A_H)
    total = 0.0
    taken = 0.0
    for A in heavy_masses
        entry = mass_yield(yields, A)
        entry === nothing && continue
        total += entry[1]
        A in on_curve && (taken += entry[1])
    end
    total > 0 || throw(
        ArgumentError(
            "the yield distribution \"$(yields.label)\" carries no positive yield over \
             $(first(heavy_masses)):$(last(heavy_masses))"
        ),
    )
    return taken / total
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
few fitted coefficients, and the method for an [`ExtractedCurve`](@ref) propagates their
covariance instead.

Throws an `ArgumentError` when the curve and the distribution share no mass number, or when the
yields summed over the shared mass numbers are not positive.
"""
function FissionFragmentsDomain.total_average(curve::RatioCurve, yields::MassYield)
    weight = Float64[]
    ratio = Float64[]
    σ_R = Float64[]
    σ_Y = Float64[]

    for (index, mass) in enumerate(curve.A_H)
        entry = mass_yield(yields, mass)
        entry === nothing && continue
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
