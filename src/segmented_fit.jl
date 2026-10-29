# Continuous piecewise-linear regression with unknown breakpoints.
#
# The model is written in the truncated-power basis
#
#   f(x) = β₀ + β₁ (x - x₀) + Σₖ γₖ (x - ψₖ)₊,   (u)₊ = max(u, 0),
#
# so continuity at every breakpoint ψₖ holds identically rather than being imposed afterwards,
# and the fit at fixed breakpoints is a single weighted linear least-squares solve which yields
# the full parameter covariance. See Muggeo, Stat. Med. 22, 3055 (2003), doi:10.1002/sim.1545, for
# the formulation.

# A design whose Gram matrix is conditioned worse than this is treated as rank-deficient.
const RANK_CONDITION_LIMIT = 1e10
# Tolerance on the unit diagonal, the symmetry and the eigenvalues of a correlation matrix.
const CORRELATION_TOLERANCE = 1e-8
# Pivots that differ by less than this, relative to their size, describe a flat curve.
const FLAT_CURVE_TOLERANCE = 1e-12

"""
    SegmentedFit

A continuous piecewise-linear fit with breakpoints selected from the data.

# Fields

- `x₀`, `x_max`: first and last abscissa of the fitted range.
- `breakpoints`: interior breakpoints `ψ`, ascending; `length(breakpoints) + 1` segments.
- `coefficients`: `[β₀,] β₁, γ₁, …` in the truncated-power basis; `β₀` is absent when the fit is
  pinned.
- `covariance`: covariance of `coefficients`: `(XᵀWX)⁻¹` for independent points and the sandwich
  `(XᵀWX)⁻¹ XᵀW Σ W X (XᵀWX)⁻¹` with `Σ = W^(-1/2) R W^(-1/2)` for points whose errors have the
  correlation matrix `R`, scaled by `max(1, wrss/expected_wrss)` where the data quote
  uncertainties and by `wrss/expected_wrss` alone where none does; see [`fit_segments`](@ref).
- `pinned_value`: the value the fit was pinned to at `x₀`, or `nothing`.
- `points`, `measured_points`: the points fitted, and how many measurements they amount to,
  `Σ f` over the measured fractions `f` of the points; equal unless the data were interpolated.
- `wrss`, `dof`, `bic`: weighted residual sum of squares `Σ (f/σ²) r²`, degrees of freedom
  `measured_points − parameters`, and the selection criterion, all counted in measurements.
  `wrss/dof` is a reduced chi-squared to the extent that the quoted uncertainties are trustworthy
  in absolute scale.
- `selection`: `(segments, bic)` for every order examined, in ascending order of segment count,
  so that the chosen order can be audited.
- `weights_imputed`: number of points that carried no uncertainty and were given the median of
  the quoted uncertainties.
- `correlated`: whether the covariance was formed with a correlation matrix of the points'
  errors, or with the points taken as independent.
- `expected_wrss`: the expectation of `wrss` the covariance scale refers to. `dof` for
  independent points; for correlated ones `dof · tr(PR)/tr(P)`, with `P` the residual projector
  of the weighted fit, smaller than `dof` for positively correlated errors because a smooth fit
  absorbs part of them.
"""
struct SegmentedFit
    x₀::Int
    x_max::Int
    breakpoints::Vector{Int}
    coefficients::Vector{Float64}
    covariance::Matrix{Float64}
    pinned_value::Union{Float64,Nothing}
    points::Int
    measured_points::Float64
    wrss::Float64
    dof::Float64
    bic::Float64
    selection::Vector{Tuple{Int,Float64}}
    weights_imputed::Int
    correlated::Bool
    expected_wrss::Float64
end

"""
    segment_count(fit) -> Int

Number of straight segments of the fit.
"""
segment_count(fit::SegmentedFit) = length(fit.breakpoints) + 1

# The truncated-power basis at `x`: `[1,] (x − x₀), (x − ψ₁)₊, …`, the intercept absent for a
# pinned fit. One row of the design matrix, and the vector a fitted value is propagated through.
function _basis_row(x::Real, x₀::Integer, ψ::AbstractVector{<:Integer}, pinned::Bool)
    row = Vector{Float64}(undef, (pinned ? 1 : 2) + length(ψ))
    index = 1
    if !pinned
        row[index] = 1.0
        index += 1
    end
    row[index] = x - x₀
    for ψₖ in ψ
        index += 1
        row[index] = max(x - ψₖ, 0)
    end
    return row
end

function _design_matrix(
    x::AbstractVector{<:Real}, x₀::Integer, ψ::AbstractVector{<:Integer}, pinned::Bool
)
    X = Matrix{Float64}(undef, length(x), (pinned ? 1 : 2) + length(ψ))
    for (i, xᵢ) in enumerate(x)
        X[i, :] .= _basis_row(xᵢ, x₀, ψ, pinned)
    end
    return X
end

"""
    fit_weights(σ) -> Tuple{Vector{Float64},Int}

Inverse-variance weights, together with the number of points whose uncertainty was not quoted.

An uncertainty is `missing` where the source quotes none. Such points carry no information about
their own weight and are given the weight of the median quoted uncertainty rather than being
discarded, so that datasets quoting no uncertainties — of which there are several — still enter
the fit; [`consensus`](@ref) imputes the same median uncertainty. This
borrows the scale of the uncertainties from the sets that do quote them, which is an assumption,
and the number of points it was applied to is recorded with the fit. Where no point has a quoted
uncertainty the weights are uniform and the fit reduces to ordinary least squares.

A quoted uncertainty must be positive. Zero denotes an exact value, which has no finite weight;
an unquoted uncertainty is `missing`, never zero, so that the two cannot be confused.

The weights are not rescaled, so that the weighted residual sum of squares is a chi-squared to
the extent that the quoted uncertainties are trustworthy. A constant factor on the weights shifts
the selection criterion by the same amount for every model order, so it does not affect which
order is chosen.

# Examples

The third point quotes no uncertainty and takes the median of the other two uncertainties,
`0.15`:

```jldoctest
julia> w, imputed = fit_weights([0.1, 0.2, missing]);

julia> (round.(w; digits = 1), imputed)
([100.0, 25.0, 44.4], 1)
```

```jldoctest
julia> fit_weights([missing, missing])
([1.0, 1.0], 2)
```
"""
function fit_weights(σ::AbstractVector{<:Union{Missing,Real}})
    quoted = .!ismissing.(σ)
    w = Vector{Float64}(undef, length(σ))
    imputed = count(!, quoted)
    if !any(quoted)
        fill!(w, 1.0)
        return (w, imputed)
    end
    for index in eachindex(σ)
        quoted[index] || continue
        σ[index] > 0 || throw(
            ArgumentError(
                "a quoted uncertainty must be positive, got $(σ[index]) at position \
                 $(index); an unquoted uncertainty is missing, and zero denotes an exact value",
            ),
        )
        w[index] = 1 / σ[index]^2
    end
    w[.!quoted] .= 1 / median(Float64[σ[i] for i in eachindex(σ) if quoted[i]])^2
    return (w, imputed)
end

# Weighted least squares at fixed breakpoints, on the weighted ordinate `yw = √w .* y`. Returns
# `nothing` where the design is rank-deficient, as when a breakpoint set leaves a segment
# without leverage.
function _solve(X::Matrix{Float64}, yw::Vector{Float64}, root_w::Vector{Float64})
    Xw = root_w .* X
    gram = Xw' * Xw
    cond(gram) < RANK_CONDITION_LIMIT || return nothing
    β = Xw \ yw
    return (β, gram)
end

# wrss over its expectation. With quoted uncertainties the weights are inverse variances, so the
# covariance is inflated where the residuals exceed the quoted scale and never shrunk below it;
# with none quoted the weights carry no scale and the ratio alone sets it.
function _covariance_scale(wrss::Float64, expected::Real, imputed::Int, n::Int)
    imputed == n && return wrss / expected
    return max(1.0, wrss / expected)
end

# A correlation matrix of the errors of n points, checked for what the covariance relies on: its
# size, its symmetry, a unit diagonal, entries within [-1, 1] and no negative eigenvalue, without
# which the variance of a fitted value could come out negative.
function _correlation_matrix(correlation::AbstractMatrix{<:Real}, n::Int)
    size(correlation) == (n, n) || throw(
        DimensionMismatch(
            "the correlation matrix of $(n) points is $(n) × $(n), got $(size(correlation))"
        ),
    )
    R = Matrix{Float64}(correlation)
    all(abs(R[m, m] - 1) ≤ CORRELATION_TOLERANCE for m in 1:n) ||
        throw(ArgumentError("the diagonal of a correlation matrix is one"))
    isapprox(R, R'; atol = CORRELATION_TOLERANCE) ||
        throw(ArgumentError("a correlation matrix is symmetric"))
    all(abs(r) ≤ 1 + CORRELATION_TOLERANCE for r in R) ||
        throw(ArgumentError("the entries of a correlation matrix lie in [-1, 1]"))
    lowest = eigmin(Symmetric(R))
    lowest ≥ -CORRELATION_TOLERANCE || throw(
        ArgumentError(
            "a correlation matrix is positive semi-definite; its smallest eigenvalue is \
             $(lowest)"
        ),
    )
    return R
end

# The covariance of the coefficients before scaling, and the expectation of wrss the scale refers
# to. For independent points these are (XᵀWX)⁻¹ and the degrees of freedom. With the weights W and
# errors of covariance Σ = W^(-1/2) R W^(-1/2), the weighted least-squares estimate keeps its
# value and its covariance is the sandwich (XᵀWX)⁻¹ XᵀWΣWX (XᵀWX)⁻¹. What a fit to such errors
# leaves of the residual sum of squares is tr(MᵀWMΣ) = tr(PR), with M = I − X(XᵀWX)⁻¹XᵀW and P
# the residual projector of the weighted design; it is referred to the fit's own degrees of
# freedom, which count measurements and breakpoints, through the ratio tr(PR)/tr(P), so that the
# two coincide for R = I.
function _coefficient_covariance(
    xs::Vector{Int}, w::Vector{Float64}, best, pinned::Bool, R::Union{Nothing,Matrix{Float64}}
)
    R === nothing && return (Symmetric(inv(best.gram)), best.dof)
    X = _design_matrix(xs, first(xs), best.ψ, pinned)
    Xw = sqrt.(w) .* X
    n = length(xs)
    p = size(X, 2)
    G⁻¹ = inv(best.gram)
    B = Xw' * R * Xw
    C = G⁻¹ * B * G⁻¹
    # R has unit diagonal, so tr(P) = n − p and tr(PR) = n − tr(G⁻¹B).
    retained = (n - tr(G⁻¹ * B)) / (n - p)
    return (Symmetric((C + C') / 2), best.dof * retained)
end

# A continuous piecewise-linear function is monotone on each segment, so its extrema over the
# range lie at the pivots: the first abscissa, the breakpoints, and the last. Checking those
# points is therefore an exact test of the bound over the whole range, not a sampled one.
function _within_bounds(
    bounds::Union{Tuple{Real,Real},Nothing},
    xs::Vector{Int},
    ψ::Vector{Int},
    β::Vector{Float64},
    pinned::Bool,
    pinned_value::Union{Real,Nothing},
)
    bounds === nothing && return true
    lower, upper = bounds
    offset = pinned_value === nothing ? 0.0 : Float64(pinned_value)
    for x in Iterators.flatten((first(xs):first(xs), ψ, last(xs):last(xs)))
        value = dot(_basis_row(x, first(xs), ψ, pinned), β) + offset
        (value > lower && value < upper) || return false
    end
    return true
end

function _bic(wrss::Float64, n::Real, parameters::Int)
    # Gaussian likelihood with the noise scale estimated from the residuals, plus the Schwarz
    # penalty, Ann. Stat. 6, 461 (1978), doi:10.1214/aos/1176344136. The breakpoints are counted
    # as parameters: they are fitted, and a criterion that ignored them would always prefer more
    # segments. On data a model fits exactly, wrss is of rounding size and the criterion
    # depends on round-off.
    return n * log(wrss / n) + parameters * log(n)
end

# Enumerate ascending interior breakpoint sets drawn from `candidates`, subject to a minimum
# number of data points per segment and, optionally, windows that must each contain a breakpoint.
#
# The first segment runs from the first abscissa to the first breakpoint and the last from the
# last breakpoint to the last abscissa, so every segment, not only the interior ones, has its
# points counted. A model without breakpoints places none in a window, so it is admissible only
# where no window is required.
function _each_breakpoint_set(
    f::Function,
    x::Vector{Int},
    candidates::Vector{Int},
    count::Int,
    min_points::Int,
    windows::Vector{UnitRange{Int}},
)
    if count == 0
        isempty(windows) && f(Int[])
        return nothing
    end
    chosen = Vector{Int}(undef, count)

    function recurse(depth::Int, start::Int)
        for index in start:length(candidates)
            ψ = candidates[index]
            lower = depth == 1 ? first(x) - 1 : chosen[depth - 1]
            Base.count(xi -> lower < xi ≤ ψ, x) ≥ min_points || continue
            chosen[depth] = ψ
            if depth == count
                Base.count(xi -> xi > ψ, x) ≥ min_points || continue
                all(window -> any(in(window), chosen), windows) || continue
                f(copy(chosen))
            else
                recurse(depth + 1, index + 1)
            end
        end
    end

    recurse(1, 1)
    return nothing
end

"""
    InsufficientDataError <: Exception

Thrown by [`fit_segments`](@ref) when the arguments are valid but the data cannot support a fit:
too few points for a single segment, fewer than three distinct abscissae, or no candidate model
that satisfies the constraints. `ArgumentError` stays reserved for invalid arguments, so a caller
can treat an unfittable dataset as an outcome without masking a programming error.
"""
struct InsufficientDataError <: Exception
    msg::String
end

function Base.showerror(io::IO, exception::InsufficientDataError)
    return print(io, "InsufficientDataError: ", exception.msg)
end

# The running best of one order of the search, mutated in place so that the closure passed to
# `_each_breakpoint_set` reassigns no captured variable.
mutable struct _BestModel
    found::Bool
    bic::Float64
    ψ::Vector{Int}
    β::Vector{Float64}
    gram::Matrix{Float64}
    wrss::Float64
end

# The best model with `count` interior breakpoints by the criterion, or `nothing` where no
# breakpoint set is admissible. Exhaustive over `candidates`.
function _best_for_order(
    xs::Vector{Int},
    ys::Vector{Float64},
    w::Vector{Float64},
    root_w::Vector{Float64},
    yw::Vector{Float64},
    candidates::Vector{Int},
    count::Int,
    min_points::Int,
    windows::Vector{UnitRange{Int}},
    pinned::Bool,
    pinned_value::Union{Real,Nothing},
    bounds::Union{Tuple{Real,Real},Nothing},
    n_measured::Float64,
)
    parameters = (pinned ? 1 : 2) + 2 * count
    dof = n_measured - parameters
    dof > 0 || return nothing
    best = _BestModel(false, NaN, Int[], Float64[], Matrix{Float64}(undef, 0, 0), NaN)
    _each_breakpoint_set(xs, candidates, count, min_points, windows) do ψ
        X = _design_matrix(xs, first(xs), ψ, pinned)
        solved = _solve(X, yw, root_w)
        solved === nothing && return nothing
        β, gram = solved
        # Counted in measurements: the weights carry each point's measured fraction.
        wrss = sum(w .* (ys .- X * β) .^ 2)
        _within_bounds(bounds, xs, ψ, β, pinned, pinned_value) || return nothing
        criterion = _bic(wrss, n_measured, parameters)
        # The first admissible candidate is kept whatever its criterion, as before.
        if !best.found || criterion < best.bic
            best.found = true
            best.bic = criterion
            best.ψ = ψ
            best.β = β
            best.gram = gram
            best.wrss = wrss
        end
        return nothing
    end
    best.found || return nothing
    return (
        ψ = best.ψ, β = best.β, gram = best.gram, wrss = best.wrss, dof = dof, bic = best.bic
    )
end

"""
    fit_segments(x, y, σ; max_segments, min_segments, min_points_per_segment, pinned_value,
                 required_windows, bounds, measured, correlation) -> SegmentedFit

Fit `y(x)` by a continuous piecewise-linear function, selecting both the number of segments and
the breakpoint positions from the data.

At a fixed set of breakpoints the fit is a weighted linear least-squares solve in the
truncated-power basis, so continuity is structural and the parameter covariance is available.
Breakpoints are restricted to abscissae that occur in the data and are searched exhaustively,
which is tractable at the size of a fragment mass range and, unlike a greedy or sequential
search, returns the global optimum of the criterion. The number of segments is chosen by the
Bayesian information criterion, which prices each additional segment and each additional
breakpoint; the criterion for every order examined is retained in the result.

The coefficient covariance is `(XᵀWX)⁻¹` for independent points and, for points whose errors
have the correlation matrix `R` given as `correlation`, the sandwich
`(XᵀWX)⁻¹ XᵀW Σ W X (XᵀWX)⁻¹` with `Σ = W^(-1/2) R W^(-1/2)`. It is scaled by `max(1, χ²/E[χ²])`
where the data quote uncertainties — inflated where the residuals exceed what the quoted
uncertainties predict, never shrunk below the quoted scale — and by `χ²/E[χ²]` alone where no
point quotes one, since uniform weights carry no scale and the noise must then be estimated from
the residuals. The expectation `E[χ²]` is `dof` for independent points and `dof · tr(PR)/tr(P)`
for correlated ones, `P` being the residual projector of the weighted fit.

# Arguments

- `max_segments`: largest number of segments examined.
- `min_segments`: smallest number examined, `1` by default. Setting it equal to `max_segments`
  fits exactly that many segments rather than selecting, which is how a given order is inspected
  on its own; the criterion is still reported for every order examined.
- `min_points_per_segment`: smallest number of data points a segment may contain; the
  identifiability guard on the search. On distinct integer abscissae a segment of `m` points
  spans at least `m − 1` units.
- `pinned_value`: when given, the fit is constrained to pass through `(first(x), pinned_value)`.
  For the multiplicity ratio of a fissioning nucleus of even mass number this is exact at the
  symmetric split, where the two fragments are identical and the ratio is one half.
- `required_windows`: mass-number windows each of which must contain a breakpoint, for imposing
  a known feature such as the minimum at the heavy magic fragment.
- `measured`: the measured fraction `f` of every point, in `(0, 1]`, all ones by default. A point
  interpolated from measurements at non-integer abscissae shares them with its neighbours and
  counts for less than one measurement: its weight is `f/σ²`, with `σ` its uncertainty as one
  measurement, and that one weight serves the solve, `χ²` and the information matrix `XᵀWX` the
  covariance is the inverse of. The sample size of the criterion and the degrees of freedom are
  `Σ f`, not the number of points. A fraction common to every point leaves the coefficients as
  they are and divides their covariance by it. Without this, correlated interpolated points
  would buy extra segments and a precision their measurements do not carry. `σ` must not carry
  the fraction already, or it enters twice.
- `correlation`: the correlation matrix `R` of the errors of the points, in their order,
  symmetric with a unit diagonal; `nothing`, the default, takes the points as independent. It
  leaves the coefficients, the breakpoints and the selected order as they are, and changes the
  covariance to the sandwich form and the reference of its scale to the expectation of `wrss`
  under `R`. An AR(1) error along the abscissa has `R[m, n] = ρ^|xₘ − xₙ|`;
  [`pooled_correlation`](@ref) gives the matrix of a curve combined from several datasets.
- `bounds`: open interval the fitted function must remain within over the whole range. For a
  ratio of the form `ν_H/(ν_L + ν_H)` the physical range is `(0, 1)`, and a fit leaving it would
  make the temperature ratio relation undefined; candidates that do are rejected outright rather
  than repaired afterwards. The test is exact, because a piecewise-linear function attains its
  extrema at its pivots.

`σ` is `missing` where no uncertainty is quoted and positive otherwise; see
[`fit_weights`](@ref).

Throws a `DimensionMismatch` when the inputs have different lengths, an `ArgumentError` for an
invalid keyword value, unsorted abscissae or a non-positive quoted uncertainty, and an
[`InsufficientDataError`](@ref) when the data cannot support a fit under the constraints.

# Examples

A ratio falling to a minimum at `A_H = 130` and rising again, recovered as two segments with the
breakpoint placed at the kink:

```jldoctest fit
julia> A_H = collect(120:139);

julia> r = [a ≤ 130 ? 0.50 - 0.020 * (a - 120) : 0.30 + 0.008 * (a - 130) for a in A_H];

julia> r .+= 0.002 .* iseven.(A_H);  # a deterministic perturbation, so the fit is not exact

julia> fit = fit_segments(A_H, r, fill(0.004, length(A_H)); max_segments = 3);

julia> segment_count(fit)
2

julia> fit.breakpoints
1-element Vector{Int64}:
 130
```

[`pivots`](@ref) returns the joined points, the form in which the fit is supplied:

```jldoctest fit
julia> [(x, round(value; digits = 3)) for (x, value) in pivots(fit)]
3-element Vector{Tuple{Int64, Float64}}:
 (120, 0.501)
 (130, 0.301)
 (139, 0.373)
```
"""
function fit_segments(
    x::AbstractVector{<:Integer},
    y::AbstractVector{<:Real},
    σ::AbstractVector{<:Union{Missing,Real}};
    max_segments::Integer = 6,
    min_segments::Integer = 1,
    min_points_per_segment::Integer = 4,
    pinned_value::Union{Real,Nothing} = nothing,
    required_windows::Vector{UnitRange{Int}} = UnitRange{Int}[],
    bounds::Union{Tuple{Real,Real},Nothing} = nothing,
    measured::Union{AbstractVector{<:Real},Nothing} = nothing,
    correlation::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
)
    length(x) == length(y) == length(σ) ||
        throw(DimensionMismatch("x, y and σ must have equal lengths, got \
                           $(length(x)), $(length(y)), $(length(σ))"))
    max_segments ≥ 1 ||
        throw(ArgumentError("max_segments must be at least 1, got $(max_segments)"))
    1 ≤ min_segments ≤ max_segments || throw(
        ArgumentError("min_segments must lie between 1 and max_segments = $(max_segments), \
             got $(min_segments)"),
    )
    min_points_per_segment ≥ 2 || throw(
        ArgumentError(
            "min_points_per_segment must be at least 2, got $(min_points_per_segment)"
        ),
    )
    issorted(x) || throw(ArgumentError("x must be sorted in ascending order"))
    R = correlation === nothing ? nothing : _correlation_matrix(correlation, length(x))

    n = length(x)
    fraction = measured === nothing ? ones(n) : collect(Float64, measured)
    length(fraction) == n || throw(
        DimensionMismatch("one measured fraction per point: $(length(fraction)) for $(n)")
    )
    all(f -> 0 < f <= 1, fraction) || throw(ArgumentError("a measured fraction lies in (0, 1]"))
    n_measured = sum(fraction)
    pinned = pinned_value !== nothing
    n ≥ min_points_per_segment + (pinned ? 0 : 1) ||
        throw(InsufficientDataError("$(n) points are too few for a single segment with \
                             min_points_per_segment = $(min_points_per_segment)"))

    xs = collect(Int, x)
    ys = pinned ? collect(Float64, y) .- Float64(pinned_value) : collect(Float64, y)
    w, imputed = fit_weights(σ)
    # The measured fraction enters here and nowhere else: one weight f/σ² for the solve, for χ²
    # and for the information matrix.
    w .*= fraction
    root_w = sqrt.(w)
    yw = root_w .* ys
    # Abscissae may repeat, as when several datasets are concatenated. Breakpoints are
    # positions, not points: the candidates are the distinct interior abscissae.
    distinct = unique(xs)
    length(distinct) ≥ 3 || throw(
        InsufficientDataError(
            "at least three distinct abscissae are needed to place a breakpoint, \
             got $(length(distinct))"
        ),
    )
    candidates = distinct[2:(end - 1)]

    best = nothing
    selection = Tuple{Int,Float64}[]

    for count in (min_segments - 1):(max_segments - 1)
        # A model with `count` interior breakpoints needs enough points for every segment, and
        # cannot use more breakpoints than there are interior abscissae.
        count ≤ length(candidates) || break
        n ≥ (count + 1) * min_points_per_segment || break

        best_for_order = _best_for_order(
            xs,
            ys,
            w,
            root_w,
            yw,
            candidates,
            Int(count),
            min_points_per_segment,
            required_windows,
            pinned,
            pinned_value,
            bounds,
            n_measured,
        )

        best_for_order === nothing && continue
        push!(selection, (count + 1, best_for_order.bic))
        if best === nothing || best_for_order.bic < best.bic
            best = best_for_order
        end
    end

    best === nothing && throw(
        InsufficientDataError(
            "no piecewise-linear model satisfies the constraints: $(n) points, \
             $(length(distinct)) distinct abscissae, \
             max_segments = $(max_segments), \
             min_points_per_segment = $(min_points_per_segment), \
             required_windows = $(required_windows), \
             pinned_value = $(pinned_value), bounds = $(bounds)"
        ),
    )

    Σ_β, expected = _coefficient_covariance(xs, w, best, pinned, R)
    scale = _covariance_scale(best.wrss, expected, imputed, n)

    return SegmentedFit(
        first(xs),
        last(xs),
        best.ψ,
        best.β,
        Matrix(Σ_β * scale),
        pinned ? Float64(pinned_value) : nothing,
        n,
        n_measured,
        best.wrss,
        best.dof,
        best.bic,
        selection,
        imputed,
        R !== nothing,
        expected,
    )
end

# The basis vector at `x`: the row of the design matrix a point there would contribute.
function _basis(fit::SegmentedFit, x::Real)
    return _basis_row(x, fit.x₀, fit.breakpoints, fit.pinned_value !== nothing)
end

"""
    evaluate(fit, x) -> Tuple{Float64,Float64}

Value and uncertainty of the fit at `x`.

The uncertainty is the propagated parameter uncertainty, `σ(x) = sqrt(jᵀ Σ j)` with `j` the basis
vector at `x`, so it widens away from the bulk of the data and vanishes at a pinned abscissa, as
it should. It is the square root of the diagonal of [`covariance`](@ref) at that abscissa; the
values at two abscissae are correlated, since they are functions of the same coefficients, and
a quantity that sums the curve over its range needs the full matrix.
"""
function evaluate(fit::SegmentedFit, x::Real)
    j = _basis(fit, x)
    value = dot(j, fit.coefficients) + something(fit.pinned_value, 0.0)
    variance = dot(j, fit.covariance, j)
    return (value, sqrt(max(variance, 0.0)))
end

"""
    covariance(fit, xs) -> Matrix{Float64}

Covariance of the fitted values at the abscissae `xs`, `J Σ Jᵀ` with `J` the rows of the basis
at each abscissa and `Σ` the coefficient covariance.

The fitted values at different abscissae are not independent — a curve with at most a handful of
coefficients cannot have independent errors at fifty mass numbers — and any quantity that combines
them, such as a yield-weighted average, propagates this matrix rather than the diagonal alone.

# Examples

The diagonal is the square of the uncertainty [`evaluate`](@ref) reports:

```jldoctest
julia> A_H = collect(120:139);

julia> r = [a ≤ 130 ? 0.50 - 0.020 * (a - 120) : 0.30 + 0.008 * (a - 130) for a in A_H];

julia> fit = fit_segments(A_H, r .+ 0.002 .* iseven.(A_H), fill(0.004, 20); max_segments = 3);

julia> Σ = covariance(fit, A_H);

julia> sqrt(Σ[11, 11]) ≈ last(evaluate(fit, 130))
true
```
"""
function covariance(fit::SegmentedFit, xs::AbstractVector{<:Real})
    J = reduce(vcat, (_basis(fit, x)' for x in xs); init = zeros(0, length(fit.coefficients)))
    Σ = J * fit.covariance * J'
    return Matrix(Symmetric(Σ))
end

"""
    evaluate(fit, x_range; label) -> RatioCurve

Evaluate the fit over a range of integer abscissae, returning a curve with propagated
uncertainties and the label `label`.
"""
function evaluate(
    fit::SegmentedFit, x_range::AbstractVector{<:Integer}; label::AbstractString = "fit"
)
    values = Float64[]
    uncertainties = Union{Missing,Float64}[]
    for x in x_range
        value, σ = evaluate(fit, Float64(x))
        push!(values, value)
        push!(uncertainties, σ)
    end
    return RatioCurve(collect(Int, x_range), values, uncertainties, String(label))
end

"""
    pivots(fit) -> Vector{Tuple{Int,Float64}}

The fit as a list of joined points, `(x, value)` at the start of the range, at every breakpoint,
and at the end.

The form a piecewise-linear curve is supplied in, and the one written to the pivot table.
"""
function pivots(fit::SegmentedFit)
    abscissae = [fit.x₀; fit.breakpoints; fit.x_max]
    return [(x, first(evaluate(fit, Float64(x)))) for x in abscissae]
end

"""
    unresolved_minimum(fit) -> Union{String,Nothing}

Why a segmented multiplicity ratio resolves no minimum, or `nothing` where it resolves one.

The multiplicity ratio falls from the symmetric split to a minimum at the heavy magic fragment and
rises beyond it. A segmented curve whose pivots show no such turn — one segment, or a range that
begins at or above the minimum — tabulates a temperature ratio without it, and the reason is
returned: a curve rising from the first mass number of its range, one falling to the last, a
flat one, or none of these. The test is of the shape of the fitted curve, an interior pivot below both of its
neighbours, and not of its `χ²`, which has no common scale across datasets.
"""
function unresolved_minimum(fit::SegmentedFit)
    points = pivots(fit)
    values = last.(points)
    m = length(values)
    any(j -> values[j] < values[j - 1] && values[j] < values[j + 1], 2:(m - 1)) &&
        return nothing
    # Flat to rounding: a fitted line through equal values has a slope of rounding size.
    maximum(values) - minimum(values) ≤ FLAT_CURVE_TOLERANCE * max(1.0, maximum(abs, values)) &&
        return "no interior minimum: the segmented ratio is flat between A_H = $(fit.x₀) and \
                $(fit.x_max)"
    lowest = argmin(values)
    lowest == 1 && return "no interior minimum: the segmented ratio rises from A_H = \
                           $(fit.x₀), the first mass number of its range"
    lowest == m && return "no interior minimum: the segmented ratio falls to A_H = \
                           $(fit.x_max), the last mass number of its range"
    return "no interior minimum between A_H = $(fit.x₀) and $(fit.x_max)"
end

"""
    unresolved_minimum(fit, windows) -> Union{String,Nothing}

Why a segmented multiplicity ratio resolves no minimum although its range covers one of the
mass-number windows a minimum lies in, or `nothing`: where it resolves one, and where its range
covers none of `windows`. A curve that begins above the window is limited in range, which its
coverage already says; one that spans the window and does not turn there tabulates a temperature
ratio without the minimum, and the reason names the window.
"""
function unresolved_minimum(fit::SegmentedFit, windows::AbstractVector{<:AbstractUnitRange})
    covered = [w for w in windows if fit.x₀ ≤ first(w) && last(w) ≤ fit.x_max]
    isempty(covered) && return nothing
    shape = unresolved_minimum(fit)
    shape === nothing && return nothing
    spans = join(("$(first(w)):$(last(w))" for w in covered), ", ")
    return "$(shape), which covers the window $(spans)"
end
