# Construction of the fragmentation range and of the isobaric charge distribution.

"""
    ChargeDistributionData

Tabulated charge polarization `ΔZ(A)` and Gaussian dispersion `rms(A)` of the isobaric charge
distribution, after Wahl, At. Data Nucl. Data Tables **38**, 1 (1988).

`fallback_ΔZ` and `fallback_rms` are used for every mass number absent from the table. The
fallback is the average behaviour quoted in the literature, `|ΔZ| = 0.5` with the sign of the
charge polarization and `rms = 0.6`; it is a coarser description than the tabulated functions,
so [`fragmentation_domain`](@ref) reports how often it was needed.
"""
struct ChargeDistributionData
    ΔZ::Dict{Int,Float64}
    rms::Dict{Int,Float64}
    fallback_ΔZ::Float64
    fallback_rms::Float64
    source::Union{String,Nothing}
end

"""
    read_charge_distribution(path, fallback_ΔZ, fallback_rms) -> ChargeDistributionData

Read a whitespace-separated table with the column layout `A ΔZ rms` and a single header line.

Passing `nothing` as `path` returns a table with no entries, so that the fallback values apply
over the whole fragmentation range. Throws an `ArgumentError` naming the file when a path is given
but cannot be read with this layout.
"""
function read_charge_distribution(
    path::Union{AbstractString,Nothing}, fallback_ΔZ::Real, fallback_rms::Real
)
    fallback_rms > 0 ||
        throw(ArgumentError("fallback rms must be positive, got $(fallback_rms)"))

    path === nothing && return ChargeDistributionData(
        Dict{Int,Float64}(), Dict{Int,Float64}(), fallback_ΔZ, fallback_rms, nothing
    )
    isfile(path) || throw(ArgumentError("charge distribution file not found: $(path)"))

    table = try
        CSV.read(
            path,
            DataFrame;
            delim = ' ',
            ignorerepeated = true,
            header = ["A", "ΔZ", "rms"],
            skipto = 2,
            types = Dict(:A => Float64, :ΔZ => Float64, :rms => Float64),
        )
    catch err
        throw(ArgumentError("charge distribution file $(path) does not have the layout \
                             `A ΔZ rms`: $(err)"))
    end

    ΔZ = Dict{Int,Float64}()
    rms = Dict{Int,Float64}()
    for row in eachrow(table)
        A = round(Int, row.A)
        ΔZ[A] = row.ΔZ
        row.rms > 0 || throw(
            ArgumentError("charge distribution file $(path) has non-positive rms at A = $(A)"),
        )
        rms[A] = row.rms
    end

    return ChargeDistributionData(ΔZ, rms, fallback_ΔZ, fallback_rms, String(path))
end

"""
    FragmentationDomain

The set of fragmentations `(A, Z)` considered, with the isobaric charge distribution `p(Z, A)`
evaluated at each of them.

For every heavy mass number `A_H` of the range, `charges_per_mass` charge numbers are taken around
the most probable charge

```
Z_p(A) = Z_UCD(A) + ΔZ(A),    Z_UCD(A) = A Z₀ / A₀,
```

and the complementary light fragment `(A₀ - A_H, Z₀ - Z_H)` is entered as well, so that the domain
covers both fragments of every pair. `p(Z, A)` is a Gaussian centred on `Z_p(A)` with dispersion
`rms(A)`, normalized within each mass number over the charges actually retained.

# Fields

- `A`, `Z`: the fragmentations, sorted by mass then charge.
- `p`: `p(Z, A)`, normalized per mass number.
- `A_H_range`: the heavy-fragment mass numbers the domain was built from.
- `fallback_masses`: mass numbers for which the tabulated `ΔZ`/`rms` were unavailable.
"""
struct FragmentationDomain
    A::Vector{Int}
    Z::Vector{Int}
    p::Vector{Float64}
    A_H_range::UnitRange{Int}
    fallback_masses::Vector{Int}
end

"""
    charges(domain, A) -> Vector{Int}

Charge numbers retained for mass number `A`, in ascending order.
"""
charges(domain::FragmentationDomain, A::Integer) = domain.Z[domain.A .== A]

"""
    charge_probability(domain, A, Z) -> Union{Float64,Missing}

The isobaric charge distribution `p(Z, A)`, or `missing` if `(A, Z)` is not in the domain.
"""
function charge_probability(domain::FragmentationDomain, A::Integer, Z::Integer)
    index = findfirst(i -> domain.A[i] == A && domain.Z[i] == Z, eachindex(domain.A))
    return index === nothing ? missing : domain.p[index]
end

"""
    most_probable_charge(A, A₀, Z₀, ΔZ) -> Float64

Most probable charge `Z_p(A) = A Z₀/A₀ + ΔZ`, the unchanged-charge-distribution value corrected
by the charge polarization.
"""
most_probable_charge(A::Integer, A₀::Integer, Z₀::Integer, ΔZ::Real) = A * Z₀ / A₀ + ΔZ

"""
    fragmentation_domain(A₀, Z₀, A_H_range, charges_per_mass, charge_data) -> FragmentationDomain

Build the fragmentation range of the fissioning nucleus `(A₀, Z₀)` over the heavy-fragment mass
numbers `A_H_range`, taking `charges_per_mass` charge numbers per mass number.

The charge polarization of a light fragment carries the opposite sign to that of its heavy
partner, which is imposed here by deriving the light fragment from its complement rather than
from the table, `Z_L = Z₀ - Z_H`.

Throws an `ArgumentError` if `charges_per_mass` is not a positive odd number — the charges are
placed symmetrically about `Z_p`, which requires an odd count — or if the range does not start at
or above the symmetric split.
"""
function fragmentation_domain(
    A₀::Integer,
    Z₀::Integer,
    A_H_range::UnitRange{Int},
    charges_per_mass::Integer,
    charge_data::ChargeDistributionData,
)
    isodd(charges_per_mass) && charges_per_mass > 0 || throw(
        ArgumentError(
            "charges_per_mass must be a positive odd number, got $(charges_per_mass)"
        ),
    )
    2 * first(A_H_range) ≥ A₀ || throw(
        ArgumentError(
            "A_H_range must start at or above the symmetric split A₀/2 = $(A₀ / 2), \
                       got $(first(A_H_range))"
        ),
    )
    last(A_H_range) < A₀ ||
        throw(ArgumentError("A_H_range must end below A₀ = $(A₀), got $(last(A_H_range))"))

    entries = Dict{Tuple{Int,Int},Float64}()
    fallback_masses = Int[]

    for A_H in A_H_range
        A_L = A₀ - A_H
        tabulated = haskey(charge_data.ΔZ, A_H) && haskey(charge_data.rms, A_H)
        tabulated || push!(fallback_masses, A_H)
        rms = tabulated ? charge_data.rms[A_H] : charge_data.fallback_rms
        # At the symmetric split the two fragments are the same nuclide, so there is nothing for
        # a charge polarization to distinguish and it must vanish. Tabulated polarizations do
        # pass through zero there; a constant fallback does not, and would centre the two charge
        # distributions of one and the same mass split on different charges, breaking the exact
        # identities R_a = 1 and R_T = 1 at symmetry. The polarization is therefore taken as zero
        # at A_H = A₀/2 whatever the source says.
        ΔZ = if 2 * A_H == A₀
            0.0
        elseif tabulated
            charge_data.ΔZ[A_H]
        else
            charge_data.fallback_ΔZ
        end

        Zₚ = most_probable_charge(A_H, A₀, Z₀, ΔZ)
        Z_H_first = round(Int, Zₚ) - (charges_per_mass - 1) ÷ 2

        for Z_H in Z_H_first:(Z_H_first + charges_per_mass - 1)
            weight = exp(-(Z_H - Zₚ)^2 / (2 * rms^2)) / (sqrt(2π) * rms)
            entries[(A_H, Z_H)] = weight
            # The symmetric split is its own complement; entering it twice would double its weight.
            if A_L != A_H
                Z_L = Z₀ - Z_H
                entries[(A_L, Z_L)] = weight
            end
        end
    end

    keys_sorted = sort!(collect(keys(entries)))
    A = [k[1] for k in keys_sorted]
    Z = [k[2] for k in keys_sorted]
    p = [entries[k] for k in keys_sorted]

    # Normalize within each mass number, over the charges actually retained.
    for mass in unique(A)
        selection = A .== mass
        total = sum(view(p, selection))
        total > 0 || throw(
            ArgumentError("isobaric charge distribution vanishes identically at A = $(mass)"),
        )
        p[selection] ./= total
    end

    return FragmentationDomain(A, Z, p, A_H_range, sort!(unique!(fallback_masses)))
end

"""
    symmetric_charge_set_is_invariant(domain, A₀, Z₀) -> Union{Bool,Missing}

Whether the charge numbers retained at the symmetric split are invariant under `Z -> Z₀ - Z`,
which is the condition for the exact identities `R_a = 1` and `R_T = 1` to hold there.

Returns `missing` when the fissioning nucleus has no symmetric split, that is, for odd `A₀`.
Invariance follows automatically from a vanishing charge polarization when `Z₀` is even; for odd
`Z₀` the most probable charge is a half-integer and rounding it to place an odd number of charges
breaks the symmetry, in which case the identities hold only approximately.
"""
function symmetric_charge_set_is_invariant(
    domain::FragmentationDomain, A₀::Integer, Z₀::Integer
)
    isodd(A₀) && return missing
    A_sym = A₀ ÷ 2
    A_sym in domain.A_H_range || return missing
    retained = Set(charges(domain, A_sym))
    return retained == Set(Z₀ - Z for Z in retained)
end

"""
    average_over_charge(values, domain, A) -> Union{Float64,Missing}

Average of `values`, a mapping from charge number to quantity, over the isobaric charge
distribution at mass number `A`.

Charges whose value is `missing` are excluded and the weights are renormalized over those that
remain, so that a quantity undefined for part of the charge distribution still yields an average
over the part where it is defined. Returns `missing` if no charge contributes.
"""
function average_over_charge(
    values::AbstractDict{Int,<:Union{Real,Missing}}, domain::FragmentationDomain, A::Integer
)
    numerator = 0.0
    denominator = 0.0
    for Z in charges(domain, A)
        value = get(values, Z, missing)
        ismissing(value) && continue
        weight = charge_probability(domain, A, Z)
        ismissing(weight) && continue
        numerator += weight * value
        denominator += weight
    end
    return denominator > 0 ? numerator / denominator : missing
end
