# The correlation, along the mass axis, of the errors of the datasets a trend is combined from.

"""
    CORRELOGRAM_LAGS

The largest lag, in mass units, of [`deviation_correlogram`](@ref).
"""
const CORRELOGRAM_LAGS = 8

# The lags an autocorrelation is fitted to where the configuration does not say.
const DEFAULT_AUTOCORRELATION_LAGS = 4

"""
    deviation_correlogram(curves, combined; lags, weights) -> Vector{Union{Missing,Float64}}

The autocorrelation, at lags `k = 1, …, lags` mass units, `lags` being [`CORRELOGRAM_LAGS`](@ref)
by default, of the deviations of the datasets `curves` from their combined curve `combined`,

```
ρₖ = Σ f d(A) d(A+k) / [Σ f d(A)² · Σ f d(A+k)²]^(1/2),
```

every sum over the pairs of mass numbers `k` apart that one dataset holds and over the datasets.
`d` is the deviation `r_ν(A_H) − r̄_ν(A_H)` of a dataset from the combined curve less its mean
over the dataset, and `f` the factor of `weights` the dataset is pooled with, one by default.

It estimates how the error of one dataset is correlated along the mass axis. A mass number held
by one dataset alone is left out: the combined value there is that dataset's own, and its
deviation vanishes whatever its error. A point without a partner `k` mass units away enters
neither sum of that lag, so a dataset on a two-unit mass grid says nothing about the odd lags and
does not lower them; what its absence at every other mass number does to the combined curve is
the business of [`pooled_correlation`](@ref). The mean of each dataset is removed, so a constant
offset of one dataset from the others is in no lag.

The estimate is low rather than high. The deviations of a short series from their own mean are
less correlated than its errors: for a series of `n` points and a lag-one correlation `ρ` the
estimate falls short by about `(1 + 3ρ)/n`. And the combined curve carries a share of every other
dataset's error into a deviation, a share that changes where one of them stops.

A dataset that shares fewer than three mass numbers with the others is left out: two deviations
less their mean are opposite and equal, a correlation of minus one whatever the errors. An entry
is `missing` at a lag no dataset spans, and every entry is `missing` where no dataset shares
three masses with another or the deviations vanish.
"""
function deviation_correlogram(
    curves::Vector{RatioCurve},
    combined::RatioCurve;
    lags::Integer = CORRELOGRAM_LAGS,
    weights::AbstractVector{<:Real} = ones(length(curves)),
)
    lags ≥ 1 || throw(ArgumentError("lags must be at least 1, got $(lags)"))
    length(weights) == length(curves) || throw(
        DimensionMismatch("one weight per curve: $(length(weights)) for $(length(curves))")
    )
    combined_at = Dict(zip(combined.A_H, combined.ratio))
    holders = Dict{Int,Int}()
    for curve in curves, A in curve.A_H
        holders[A] = get(holders, A, 0) + 1
    end
    shared(A) = haskey(combined_at, A) && holders[A] ≥ 2

    lagged = zeros(lags)
    leading = zeros(lags)
    trailing = zeros(lags)
    squared = 0.0
    terms = 0
    for (curve, f) in zip(curves, weights)
        masses = [A for A in curve.A_H if shared(A)]
        # Two points less their mean are opposite and equal: a lag of minus one by construction.
        length(masses) < 3 && continue
        deviation = Dict(
            A => r - combined_at[A] for (A, r) in zip(curve.A_H, curve.ratio) if shared(A)
        )
        offset = mean(values(deviation))
        for A in masses
            d = deviation[A] - offset
            terms += 1
            squared += d^2
            for k in 1:lags
                haskey(deviation, A + k) || continue
                partner = deviation[A + k] - offset
                lagged[k] += f * d * partner
                leading[k] += f * d^2
                trailing[k] += f * partner^2
            end
        end
    end
    # Deviations at the level of rounding, as between identical datasets, say nothing.
    correlogram = Vector{Union{Missing,Float64}}(missing, lags)
    (terms == 0 || sqrt(squared / terms) ≤ 1e-12) && return correlogram
    for k in 1:lags
        scale = sqrt(leading[k] * trailing[k])
        scale > 0 && (correlogram[k] = lagged[k] / scale)
    end
    return correlogram
end

"""
    autocorrelation_decay(correlogram; lags = 4) -> Float64

The coefficient `ρ` of the AR(1) kernel `ρ^k` that follows the first `lags` lags of a correlogram
in least squares: the `ρ` minimising `Σₖ (ρₖ − ρ^k)²` over the lags `k ≤ lags` the correlogram
holds a value at, `0.0` where it holds none.

It is fitted to several lags rather than read from the first: the decay of the measured
correlograms is slower than the powers of their first lag. It is bounded below by zero: a
negative estimate is taken as no correlation, since an anticorrelation would bring the covariance
of a smooth fit below that of independent points on the strength of a noisy estimate. It is
bounded above by 0.999, short of the kernel of a fully correlated error.

The minimum is located on a grid of step 0.001 over `[0, 0.999]` and refined by golden-section
search between the neighbours of the grid minimum, so the result is deterministic.
"""
function autocorrelation_decay(
    correlogram::AbstractVector{<:Union{Missing,Real}};
    lags::Integer = DEFAULT_AUTOCORRELATION_LAGS,
)
    lags ≥ 1 || throw(ArgumentError("lags must be at least 1, got $(lags)"))
    orders = Int[]
    estimates = Float64[]
    for k in 1:min(lags, length(correlogram))
        ρₖ = correlogram[k]
        ismissing(ρₖ) && continue
        push!(orders, k)
        push!(estimates, ρₖ)
    end
    isempty(orders) && return 0.0
    S(ρ) = sum((ρₖ - ρ^k)^2 for (k, ρₖ) in zip(orders, estimates))

    grid = 0.0:0.001:0.999
    i = argmin(S.(grid))
    i == 1 && return 0.0
    a = clamp(grid[max(i - 1, 1)], 0.0, 0.999)
    b = clamp(grid[min(i + 1, length(grid))], 0.0, 0.999)
    φ = (sqrt(5) - 1) / 2
    c = b - φ * (b - a)
    d = a + φ * (b - a)
    S_c = S(c)
    S_d = S(d)
    for _ in 1:60
        if S_c ≤ S_d
            b, d, S_d = d, c, S_c
            c = b - φ * (b - a)
            S_c = S(c)
        else
            a, c, S_c = c, d, S_d
            d = a + φ * (b - a)
            S_d = S(d)
        end
    end
    return clamp((a + b) / 2, 0.0, 0.999)
end

"""
    pooled_correlation(A_H, members, shares, ρ) -> Matrix{Float64}

The correlation matrix of the errors of a combined curve whose datasets each err along the mass
axis with an AR(1) correlation `ρ` and independently of one another:

```
R(A, A′) = ρ^|A − A′| Σᵢ [sᵢ(A) sᵢ(A′)]^(1/2),
```

the sum over the datasets `i` that hold both mass numbers, `sᵢ(A)` being the share of the
combined value at `A` that dataset `i` carries, its normalized weight. `members[m]` lists the
datasets combined at `A_H[m]` and `shares[m]` their shares, which sum to one.

Where the same datasets hold two mass numbers with the same shares the entry is `ρ^|A − A′|`. A
dataset held at one of them and not at the other brings an error the other combined value does
not share, and the entry is lower by its part: a dataset on a two-unit mass grid lowers the
correlation between neighbouring combined points and not that between next-neighbours. The
diagonal is one, and the matrix is positive semi-definite, being the element-wise product of two
such matrices.
"""
function pooled_correlation(
    A_H::AbstractVector{<:Integer},
    members::AbstractVector{<:AbstractVector{<:Integer}},
    shares::AbstractVector{<:AbstractVector{<:Real}},
    ρ::Real,
)
    n = length(A_H)
    length(members) == n == length(shares) ||
        throw(DimensionMismatch("one list of members and one of shares per mass number: \
             $(length(members)) and $(length(shares)) for $(n)"))
    0 ≤ ρ < 1 || throw(ArgumentError("ρ must lie in [0, 1), got $(ρ)"))
    for m in 1:n
        length(members[m]) == length(shares[m]) ||
            throw(DimensionMismatch("one share per member at A_H = $(A_H[m])"))
        isapprox(sum(shares[m]), 1; atol = 1e-8) || throw(
            ArgumentError("the shares at A_H = $(A_H[m]) sum to $(sum(shares[m])), not to one"),
        )
    end
    share_at = [Dict{Int,Float64}(zip(members[m], shares[m])) for m in 1:n]
    R = Matrix{Float64}(undef, n, n)
    for m in 1:n
        R[m, m] = 1.0
        for k in (m + 1):n
            overlap = 0.0
            for (i, s) in share_at[m]
                overlap += sqrt(s * get(share_at[k], i, 0.0))
            end
            R[m, k] = R[k, m] = ρ^abs(A_H[m] - A_H[k]) * overlap
        end
    end
    return R
end
