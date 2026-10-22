abstract type ADSEstimator end

"""
    IntervalBound()

Use outward-rounded Taylor-model remainders and retained coefficient widths as
uniform absolute error bounds in `adaptive_map`. Requires IntervalArithmetic.
The original function is reevaluated on every child. Supports static maps,
positive absolute tolerances and no guard degrees or sampled acceptance checks.
The returned map retains its interval coefficients, remainder and domain.
Domain shape and split direction are independent of the error method.
Analytic domain failures trigger subdivision. A failure remaining at a splitting
limit throws even with `strict=false`, since no valid model exists for that cell.
"""
struct IntervalBound <: ADSEstimator end

"""
    GuardedTail()

Estimate truncation error using extra Taylor degrees in [`adaptive_map`](@ref)
and [`adaptive_flow`](@ref). Sum the absolute discarded coefficients and
extrapolate their degree norms to the next degree. The default `guard_order`
is two. An exactly represented polynomial has zero estimated error.
Point checks supplement this estimate by default. It is not a rigorous bound
on the uncomputed Taylor remainder.
"""
struct GuardedTail <: ADSEstimator end

"""
    ExtrapolatedTail()

Fit `log(maximum absolute coefficient)` against total degree, omitting degrees
below two and zero coefficients, and extrapolate to `order + 1`. With fewer than
two nonzero nonlinear degrees, use the largest nonlinear coefficient in the last
two degrees. This ADS indicator requires no guard degrees, but may over-refine
exact polynomials or miss terms absent from the computed expansion.
At order one the estimate is zero;
point checks, enabled by default, provide the acceptance check.
"""
struct ExtrapolatedTail <: ADSEstimator end

"""
    LastTerms(; degrees = 2, safety = 1)

Use `safety` times the largest absolute coefficient in the last `degrees`
retained total degrees above the affine part as an ADS error indicator. Two degrees
accommodate even/odd series. `degrees` must be positive and `safety` finite and at least
one. Coefficients already include the physical box radii through normalized
coordinates; they are not scaled a second time. This inexpensive indicator
requires no guard degrees and is not a remainder bound. At order one its estimate
is zero, so use the default point checks to check acceptance.
"""
struct LastTerms{T <: Real} <: ADSEstimator
    degrees::Int
    safety::T
    function LastTerms(; degrees::Integer = 2, safety::Real = 1)
        1 <= degrees <= typemax(Int) || throw(ArgumentError("degrees must be positive"))
        isfinite(safety) && safety >= 1 || throw(ArgumentError("safety must be finite and at least one"))
        return new{typeof(safety)}(Int(degrees), safety)
    end
end

# Fit log(norm) against degree without a matrix or temporary logarithm array.
# The target is one degree beyond the last element of norms.
function ads_next_norm(norms; fallback = zero(eltype(norms)), first_degree = 1)
    R = eltype(norms)
    count = 0; sx = sy = sxx = sxy = zero(R)
    for d in first_degree:(length(norms) - 1)
        c = norms[d + 1]
        iszero(c) && continue
        count += 1
        y = log(c)
        sx += d; sy += y; sxx += R(d)^2; sxy += d * y
    end
    count < 2 && return fallback
    slope = (count * sxy - sx * sy) / (count * sxx - sx * sx)
    return exp((sy - slope * sx) / count + slope * length(norms))
end

function ads_error(p, order, ::GuardedTail)
    norms = @view degree_norms(p, 0, 1)[1:(working_order(p.algebra) + 1)]
    tail = sum(@view norms[(order + 2):end])
    return iszero(tail) ? tail : tail + ads_next_norm(norms)
end
function ads_error(p, order, ::ExtrapolatedTail)
    norms = degree_norms(p, 0, 0)
    retained = @view norms[1:(order + 1)]
    fallback = maximum(
        @view(retained[(max(2, order - 1) + 1):end]);
        init = zero(eltype(norms))
    )
    return ads_next_norm(retained; fallback, first_degree = 2)
end
function ads_error(p, order, estimator::LastTerms)
    norms = degree_norms(p, 0, 0)
    tail = maximum(
        @view(norms[(max(2, order - estimator.degrees + 1) + 1):(order + 1)]);
        init = zero(eltype(norms))
    )
    safety = oftype(tail, estimator.safety)
    isfinite(safety) || throw(ArgumentError("safety is not representable in the coefficient type"))
    return safety * tail
end

ads_tail_start(order, ::GuardedTail) = order + 1
ads_tail_start(order, ::ExtrapolatedTail) = max(2, order - 1)
ads_tail_start(order, estimator::LastTerms) = max(2, order - estimator.degrees + 1)

# Expected reduction in the selected tail when normalized coordinate i is halved.
function ads_contributions!(scores, p, row, order, estimator)
    basis = valid(p).basis
    C = eltype(scores)
    first = basis.ends[ads_tail_start(order, estimator)] + 1
    for k in first:p.len
        c = abs(p.coeffs[k])
        iszero(c) && continue
        for i in axes(scores, 2)
            scores[row, i] += c * (one(C) - exp2(-C(basis.exponents[i, k])))
        end
    end
    return scores
end
