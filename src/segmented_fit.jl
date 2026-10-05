# Continuous piecewise-linear regression with unknown breakpoints.
#
# The multiplicity ratio r_ν(A_H) has a known systematic shape — below one half from the symmetric
# split to the most probable fragmentation, a minimum at the heavy magic fragment, then an almost
# linear rise — while the ratio extracted point by point from ν(A) data is scattered and, for some
# datasets, sparse. Describing the ratio by a small number of joined straight segments imposes
# that shape, and the temperature ratio obtained from the fitted ratio is smooth enough to serve
# as tabulated input elsewhere.
#
# The model is written in the truncated-power basis
#
#   f(x) = β₀ + β₁ (x - x₀) + Σₖ γₖ (x - ψₖ)₊,   (u)₊ = max(u, 0),
#
# so continuity at every breakpoint ψₖ holds identically rather than being imposed afterwards,
# and the fit at fixed breakpoints is a single weighted linear least-squares solve which yields
# the full parameter covariance. See Muggeo, Stat. Med. 22, 3055 (2003), doi:10.1002/sim.1545, for
# the formulation.

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
- `weights_imputed`: number of points that carried no uncertainty and were given the median
  weight.
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
    segments(fit) -> Int

Number of straight segments of the fit.
"""
segments(fit::SegmentedFit) = length(fit.breakpoints) + 1

function _design_matrix(
    x::AbstractVector{<:Real}, x₀::Integer, ψ::AbstractVector{<:Integer}, pinned::Bool
)
    n = length(x)
    columns = Vector{Vector{Float64}}()
    pinned || push!(columns, ones(n))
    push!(columns, Float64[xi - x₀ for xi in x])
    for ψₖ in ψ
        push!(columns, Float64[max(xi - ψₖ, 0) for xi in x])
    end
    return reduce(hcat, columns)
end

"""
    fit_weights(σ) -> Tuple{Vector{Float64},Int}

Inverse-variance weights, together with the number of points whose uncertainty was not quoted.

An uncertainty is `missing` where the source quotes none. Such points carry no information about
their own weight and are given the median of the quoted weights rather than being discarded, so
that datasets quoting no uncertainties — of which there are several — still enter the fit. This
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

The third point quotes no uncertainty and takes the median of the other two weights:

```jldoctest
julia> w, imputed = fit_weights([0.1, 0.2, missing]);

julia> (round.(w; digits = 1), imputed)
([100.0, 25.0, 62.5], 1)
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
    w[.!quoted] .= median(view(w, quoted))
    return (w, imputed)
end

# Weighted least squares at fixed breakpoints. Returns nothing if the design is rank-deficient,
# which happens when a candidate breakpoint set leaves a segment without enough leverage.
function _solve(X::Matrix{Float64}, y::Vector{Float64}, w::Vector{Float64})
    root_w = sqrt.(w)
    Xw = root_w .* X
    yw = root_w .* y
    gram = Xw' * Xw
    (isfinite(cond(gram)) && cond(gram) < 1e10) || return nothing
    β = Xw \ yw
    residual = yw .- Xw * β
    wrss = sum(abs2, residual)
    return (β, wrss, gram)
end

# The factor the coefficient covariance is scaled by: wrss over its expectation, which is the
# degrees of freedom for independent points. Where the data quote uncertainties the weights are
# inverse variances and the ratio is a reduced chi-squared: the covariance is inflated where the
# residuals exceed what the quoted uncertainties predict and left alone where they do not, since
# shrinking it below the quoted scale would claim a precision the data do not carry. Where no
# point quotes an uncertainty the weights are uniform and carry no scale at all; the noise is then
# estimated from the residuals, which is ordinary least squares.
function _covariance_scale(wrss::Float64, expected::Real, imputed::Int, n::Int)
    imputed == n && return wrss / expected
    return max(1.0, wrss / expected)
end

# A correlation matrix of the errors of n points, checked for what the covariance relies on: its
# size, its symmetry, a unit diagonal and entries within [-1, 1].
function _correlation_matrix(correlation::AbstractMatrix{<:Real}, n::Int)
    size(correlation) == (n, n) || throw(
        DimensionMismatch(
            "the correlation matrix of $(n) points is $(n) × $(n), got $(size(correlation))"
        ),
    )
    R = Matrix{Float64}(correlation)
    tolerance = 1e-8
    all(abs(R[m, m] - 1) ≤ tolerance for m in 1:n) ||
        throw(ArgumentError("the diagonal of a correlation matrix is one"))
    isapprox(R, R'; atol = tolerance) ||
        throw(ArgumentError("a correlation matrix is symmetric"))
    all(abs(r) ≤ 1 + tolerance for r in R) ||
        throw(ArgumentError("the entries of a correlation matrix lie in [-1, 1]"))
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
        value = offset
        index = 1
        if !pinned
            value += β[index]
            index += 1
        end
        value += β[index] * (x - first(xs))
        index += 1
        for (k, ψₖ) in enumerate(ψ)
            value += β[index + k - 1] * max(x - ψₖ, 0)
        end
        (value > lower && value < upper) || return false
    end
    return true
end

function _bic(wrss::Float64, n::Real, parameters::Int)
    # Gaussian likelihood with the noise scale estimated from the residuals, plus the Schwarz
    # penalty, Ann. Stat. 6, 461 (1978), doi:10.1214/aos/1176344136. The breakpoints are counted
    # as parameters: they are fitted, and a criterion that ignored them would always prefer more
    # segments.
    return n * log(wrss / n) + parameters * log(n)
end

# Enumerate ascending interior breakpoint sets drawn from `candidates`, subject to a minimum
# number of data points per segment, a minimum extent per segment in units of the abscissa, and,
# optionally, windows that must each contain a breakpoint.
#
# The first segment runs from the first abscissa to the first breakpoint and the last from the
# last breakpoint to the last abscissa, so every segment has an extent to test, not only the
# interior ones. A model without breakpoints places none in a window, so it is admissible only
# where no window is required.
function _each_breakpoint_set(
    f::Function,
    x::Vector{Int},
    candidates::Vector{Int},
    count::Int,
    min_points::Int,
    min_span::Int,
    windows::Vector{UnitRange{Int}},
)
    if count == 0
        isempty(windows) && last(x) - first(x) ≥ min_span && f(Int[])
        return nothing
    end
    chosen = Vector{Int}(undef, count)

    function recurse(depth::Int, start::Int)
        for index in start:length(candidates)
            ψ = candidates[index]
            lower = depth == 1 ? first(x) - 1 : chosen[depth - 1]
            Base.count(xi -> lower < xi ≤ ψ, x) ≥ min_points || continue
            ψ - (depth == 1 ? first(x) : chosen[depth - 1]) ≥ min_span || continue
            chosen[depth] = ψ
            if depth == count
                Base.count(xi -> xi > ψ, x) ≥ min_points || continue
                last(x) - ψ ≥ min_span || continue
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

"""
    fit_segments(x, y, σ; max_segments, min_points_per_segment, min_segment_span, pinned_value,
                 required_windows, bounds) -> SegmentedFit

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
  identifiability guard on the search.
- `min_segment_span`: smallest extent of a segment in units of the abscissa, from its first
  pivot to its last, `1` by default. The point count bounds how much data a segment rests on and
  the span bounds how short a feature it may assert; at four points per segment on consecutive
  mass numbers the two coincide at three mass units, so the span guard acts where the abscissae
  repeat or the point count is set lower.
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

julia> segments(fit)
2

julia> fit.breakpoints
1-element Vector{Int64}:
 130
```

[`pivots`](@ref) returns the joined points, which is the form the parameterization is supplied in:

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
    min_segment_span::Integer = 1,
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
    min_segment_span ≥ 1 ||
        throw(ArgumentError("min_segment_span must be at least 1, got $(min_segment_span)"))
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
    ys = collect(Float64, y)
    pinned && (ys = ys .- pinned_value)
    w, imputed = fit_weights(σ)
    # The measured fraction enters here and nowhere else: one weight f/σ² for the solve, for χ²
    # and for the information matrix.
    w .*= fraction
    # Several datasets of one fissioning nucleus contribute points at the same mass number, so
    # the abscissae repeat. Breakpoints are positions, not points: the candidates are the distinct
    # interior mass numbers.
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

        best_for_order = nothing
        _each_breakpoint_set(
            xs,
            candidates,
            count,
            min_points_per_segment,
            Int(min_segment_span),
            required_windows,
        ) do ψ
            X = _design_matrix(xs, first(xs), ψ, pinned)
            solved = _solve(X, ys, w)
            solved === nothing && return nothing
            β, _, gram = solved
            parameters = length(β) + length(ψ)
            # Counted in measurements: the weights carry each point's measured fraction.
            wrss = sum(w .* (ys .- X * β) .^ 2)
            dof = n_measured - parameters
            dof > 0 || return nothing
            _within_bounds(bounds, xs, ψ, β, pinned, pinned_value) || return nothing
            criterion = _bic(wrss, n_measured, parameters)
            if best_for_order === nothing || criterion < best_for_order.bic
                best_for_order = (
                    ψ = ψ, β = β, gram = gram, wrss = wrss, dof = dof, bic = criterion
                )
            end
            return nothing
        end

        best_for_order === nothing && continue
        push!(selection, (count + 1, best_for_order.bic))
        if best === nothing || best_for_order.bic < best.bic
            best = best_for_order
        end
    end

    best === nothing && throw(
        InsufficientDataError(
            "no piecewise-linear model satisfies the constraints: $(n) points, \
             max_segments = $(max_segments), \
             min_points_per_segment = $(min_points_per_segment), \
             min_segment_span = $(min_segment_span), \
             required_windows = $(required_windows)"
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
    j = Float64[]
    fit.pinned_value === nothing && push!(j, 1.0)
    push!(j, x - fit.x₀)
    for ψ in fit.breakpoints
        push!(j, max(x - ψ, 0))
    end
    return j
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
    evaluate(fit, x_range) -> RatioCurve

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

This is the representation a piecewise-linear parameterization is normally supplied in, and the
form written to the segment output file.
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
returned: a curve rising from the first mass number of its range, one falling to the last, or
neither. The test is of the shape of the fitted curve, an interior pivot below both of its
neighbours, and not of its `χ²`, which has no common scale across datasets.
"""
function unresolved_minimum(fit::SegmentedFit)
    points = pivots(fit)
    values = last.(points)
    m = length(values)
    any(j -> values[j] < values[j - 1] && values[j] < values[j + 1], 2:(m - 1)) &&
        return nothing
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
