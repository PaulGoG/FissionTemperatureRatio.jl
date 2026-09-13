# Continuous piecewise-linear regression with unknown breakpoints.
#
# The multiplicity ratio r_ν(A_H) has a known systematic shape — below one half from the symmetric
# split to the most probable fragmentation, a minimum at the heavy magic fragment, then an almost
# linear rise — while the ratio extracted point by point from ν(A) data is scattered and, for some
# data sets, sparse. Describing the ratio by a small number of joined straight segments imposes
# that shape, and the temperature ratio obtained from the fitted ratio is smooth enough to serve
# as tabulated input elsewhere.
#
# The model is written in the truncated-power basis
#
#   f(x) = β₀ + β₁ (x - x₀) + Σₖ γₖ (x - ψₖ)₊,   (u)₊ = max(u, 0),
#
# so continuity at every breakpoint ψₖ holds identically rather than being imposed afterwards,
# and the fit at fixed breakpoints is a single weighted linear least-squares solve which yields
# the full parameter covariance. See Muggeo, Stat. Med. 22, 3055 (2003) for the formulation.

"""
    SegmentedFit

A continuous piecewise-linear fit with breakpoints selected from the data.

# Fields

- `x₀`, `x_max`: first and last abscissa of the fitted range.
- `breakpoints`: interior breakpoints `ψ`, ascending; `length(breakpoints) + 1` segments.
- `coefficients`: `[β₀,] β₁, γ₁, …` in the truncated-power basis; `β₀` is absent when the fit is
  pinned.
- `covariance`: covariance of `coefficients`, scaled by the reduced chi-squared.
- `pinned_value`: the value the fit was pinned to at `x₀`, or `nothing`.
- `wrss`, `dof`, `bic`: weighted residual sum of squares, degrees of freedom, and the
  selection criterion. `wrss/dof` is a reduced chi-squared to the extent that the quoted
  uncertainties are trustworthy in absolute scale.
- `selection`: `(segments, bic)` for every order examined, in ascending order of segment count,
  so that the chosen order can be audited.
- `weights_imputed`: number of points that carried no uncertainty and were given the median
  weight.
"""
struct SegmentedFit
    x₀::Int
    x_max::Int
    breakpoints::Vector{Int}
    coefficients::Vector{Float64}
    covariance::Matrix{Float64}
    pinned_value::Union{Float64,Nothing}
    wrss::Float64
    dof::Int
    bic::Float64
    selection::Vector{Tuple{Int,Float64}}
    weights_imputed::Int
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

Inverse-variance weights, together with the number of points whose uncertainty was absent.

Points without an uncertainty carry no information about their own weight. They are given the
median of the positive weights rather than being discarded, so that data sets quoting no
uncertainties — of which there are several — still enter the fit. This borrows the scale of the
uncertainties from the sets that do quote them, which is an assumption, and the number of points
it was applied to is recorded with the fit. Where no point has an uncertainty the weights are
uniform and the fit reduces to ordinary least squares.

The weights are not rescaled, so that the weighted residual sum of squares is a chi-squared to
the extent that the quoted uncertainties are trustworthy. A constant factor on the weights shifts
the selection criterion by the same amount for every model order, so it does not affect which
order is chosen.

# Examples

The third point quotes no uncertainty and takes the median of the other two weights:

```jldoctest
julia> w, imputed = fit_weights([0.1, 0.2, 0.0]);

julia> (round.(w; digits = 1), imputed)
([100.0, 25.0, 62.5], 1)
```

```jldoctest
julia> fit_weights([0.0, 0.0])
([1.0, 1.0], 2)
```
"""
function fit_weights(σ::AbstractVector{<:Real})
    positive = σ .> 0
    w = Vector{Float64}(undef, length(σ))
    imputed = count(!, positive)
    if !any(positive)
        fill!(w, 1.0)
        return (w, imputed)
    end
    w[positive] .= 1 ./ σ[positive] .^ 2
    w[.!positive] .= median(view(w, positive))
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
    offset = pinned ? Float64(pinned_value) : 0.0
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

function _bic(wrss::Float64, n::Int, parameters::Int)
    # Gaussian likelihood with the noise scale estimated from the residuals, plus the Schwarz
    # penalty, Ann. Stat. 6, 461 (1978). The breakpoints are counted as parameters: they are
    # fitted, and a criterion that ignored them would always prefer more segments.
    return n * log(wrss / n) + parameters * log(n)
end

# Enumerate ascending interior breakpoint sets drawn from `candidates`, subject to a minimum
# number of data points per segment and, optionally, windows that must each contain a breakpoint.
function _each_breakpoint_set(
    f::Function,
    x::Vector{Int},
    candidates::Vector{Int},
    count::Int,
    min_points::Int,
    windows::Vector{UnitRange{Int}},
)
    count == 0 && return f(Int[])
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
    fit_segments(x, y, σ; max_segments, min_points_per_segment, pinned_value,
                 required_windows) -> SegmentedFit

Fit `y(x)` by a continuous piecewise-linear function, selecting both the number of segments and
the breakpoint positions from the data.

At a fixed set of breakpoints the fit is a weighted linear least-squares solve in the
truncated-power basis, so continuity is structural and the parameter covariance is available.
Breakpoints are restricted to abscissae that occur in the data and are searched exhaustively,
which is tractable at the size of a fragment mass range and, unlike a greedy or sequential
search, returns the global optimum of the criterion. The number of segments is chosen by the
Bayesian information criterion, which prices each additional segment and each additional
breakpoint; the criterion for every order examined is retained in the result.

# Arguments

- `max_segments`: largest number of segments examined.
- `min_points_per_segment`: smallest number of data points a segment may contain; the
  identifiability guard on the search.
- `pinned_value`: when given, the fit is constrained to pass through `(first(x), pinned_value)`.
  For the multiplicity ratio of a fissioning nucleus of even mass number this is exact at the
  symmetric split, where the two fragments are identical and the ratio is one half.
- `required_windows`: mass-number windows each of which must contain a breakpoint, for imposing
  a known feature such as the minimum at the heavy magic fragment.
- `bounds`: open interval the fitted function must remain within over the whole range. For a
  ratio of the form `ν_H/(ν_L + ν_H)` the physical range is `(0, 1)`, and a fit leaving it would
  make the temperature ratio relation undefined; candidates that do are rejected outright rather
  than repaired afterwards. The test is exact, because a piecewise-linear function attains its
  extrema at its pivots.

Throws an `ArgumentError` when the inputs have different lengths, when fewer points are available
than the smallest model requires, or when no candidate model satisfies the constraints.

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
    σ::AbstractVector{<:Real};
    max_segments::Integer = 6,
    min_points_per_segment::Integer = 4,
    pinned_value::Union{Real,Nothing} = nothing,
    required_windows::Vector{UnitRange{Int}} = UnitRange{Int}[],
    bounds::Union{Tuple{Real,Real},Nothing} = nothing,
)
    length(x) == length(y) == length(σ) ||
        throw(DimensionMismatch("x, y and σ must have equal lengths, got \
                           $(length(x)), $(length(y)), $(length(σ))"))
    max_segments ≥ 1 ||
        throw(ArgumentError("max_segments must be at least 1, got $(max_segments)"))
    min_points_per_segment ≥ 2 || throw(
        ArgumentError(
            "min_points_per_segment must be at least 2, got $(min_points_per_segment)"
        ),
    )
    issorted(x) || throw(ArgumentError("x must be sorted in ascending order"))

    n = length(x)
    pinned = pinned_value !== nothing
    n ≥ min_points_per_segment + (pinned ? 0 : 1) ||
        throw(ArgumentError("$(n) points are too few for a single segment with \
                             min_points_per_segment = $(min_points_per_segment)"))

    xs = collect(Int, x)
    ys = collect(Float64, y)
    pinned && (ys = ys .- pinned_value)
    w, imputed = fit_weights(σ)
    # Several data sets of one fissioning nucleus contribute points at the same mass number, so
    # the abscissae repeat. Breakpoints are positions, not points: the candidates are the distinct
    # interior mass numbers.
    distinct = unique(xs)
    length(distinct) ≥ 3 || throw(
        ArgumentError("at least three distinct abscissae are needed to place a breakpoint, \
                       got $(length(distinct))"),
    )
    candidates = distinct[2:(end - 1)]

    best = nothing
    selection = Tuple{Int,Float64}[]

    for count in 0:(max_segments - 1)
        # A model with `count` interior breakpoints needs enough points for every segment, and
        # cannot use more breakpoints than there are interior abscissae.
        count ≤ length(candidates) || break
        n ≥ (count + 1) * min_points_per_segment || break

        best_for_order = nothing
        _each_breakpoint_set(
            xs, candidates, count, min_points_per_segment, required_windows
        ) do ψ
            X = _design_matrix(xs, first(xs), ψ, pinned)
            solved = _solve(X, ys, w)
            solved === nothing && return nothing
            β, wrss, gram = solved
            parameters = length(β) + length(ψ)
            dof = n - parameters
            dof > 0 || return nothing
            _within_bounds(bounds, xs, ψ, β, pinned, pinned_value) || return nothing
            criterion = _bic(wrss, n, parameters)
            if best_for_order === nothing || criterion < best_for_order.bic
                covariance = Symmetric(inv(gram)) * (wrss / dof)
                best_for_order = (
                    ψ = ψ,
                    β = β,
                    covariance = Matrix(covariance),
                    wrss = wrss,
                    dof = dof,
                    bic = criterion,
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
        ArgumentError("no piecewise-linear model satisfies the constraints: $(n) points, \
                       max_segments = $(max_segments), \
                       min_points_per_segment = $(min_points_per_segment), \
                       required_windows = $(required_windows)")
    )

    return SegmentedFit(
        first(xs),
        last(xs),
        best.ψ,
        best.β,
        best.covariance,
        pinned ? Float64(pinned_value) : nothing,
        best.wrss,
        best.dof,
        best.bic,
        selection,
        imputed,
    )
end

"""
    evaluate(fit, x) -> Tuple{Float64,Float64}

Value and uncertainty of the fit at `x`.

The uncertainty is the propagated parameter uncertainty, `σ(x) = sqrt(jᵀ Σ j)` with `j` the basis
vector at `x`, so it widens away from the bulk of the data and vanishes at a pinned abscissa, as
it should. The covariance is scaled by the reduced chi-squared, which estimates the noise level
from the residuals rather than trusting the quoted uncertainties in absolute terms.
"""
function evaluate(fit::SegmentedFit, x::Real)
    pinned = fit.pinned_value !== nothing
    j = Float64[]
    pinned || push!(j, 1.0)
    push!(j, x - fit.x₀)
    for ψ in fit.breakpoints
        push!(j, max(x - ψ, 0))
    end
    value = dot(j, fit.coefficients) + (pinned ? fit.pinned_value : 0.0)
    variance = dot(j, fit.covariance, j)
    return (value, sqrt(max(variance, 0.0)))
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
    uncertainties = Float64[]
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
