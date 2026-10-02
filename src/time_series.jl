"""
    TimeSeries(coefficients::AbstractVector{<:Real})

A truncated expansion `sum(coefficients[k+1] * τ^k)` in a local time offset.
Coefficients may themselves be [`TaylorPolynomial`](@ref)s: the time order
and the polynomial algebra's uncertainty order are independent. Construction
copies the coefficients. No global algebra is needed for numeric coefficients.

Use `coefficient(s, k)` (starting at zero), `max_order(s)`, `s(τ)`,
`differentiate(s)` and `integrate(s)` to inspect or evaluate the series.
Arithmetic between series requires equal orders; scalar arithmetic preserves
the order. Elementary powers, roots, exponential, logarithmic, trigonometric
and hyperbolic functions are supported at regular expansion points.
"""
struct TimeSeries{T <: Real} <: Number
    coefficients::Vector{T}
    function TimeSeries{T}(coefficients::Vector{T}) where {T <: Real}
        isempty(coefficients) && throw(ArgumentError("A time series needs a constant coefficient"))
        return new{T}(coefficients)
    end
end
TimeSeries(c::AbstractVector{T}) where {T <: Real} = TimeSeries{T}(copy.(c))
coefficient_type(::TimeSeries{T}) where {T} = T
max_order(s::TimeSeries) = length(s.coefficients) - 1
constant_term(s::TimeSeries) = first(s.coefficients)
function coefficient(s::TimeSeries, k::Integer)
    0 <= k <= max_order(s) || throw(BoundsError(s.coefficients, k + 1))
    return s.coefficients[k + 1]
end
Base.copy(s::TimeSeries) = TimeSeries(s.coefficients)
Base.iszero(s::TimeSeries) = all(iszero, s.coefficients)
Base.isfinite(s::TimeSeries) = all(isfinite, s.coefficients)
Base.real(s::TimeSeries) = s
Base.conj(s::TimeSeries) = s
Base.abs2(s::TimeSeries) = s * s
Base.:+(s::TimeSeries) = s
Base.:-(s::TimeSeries) = TimeSeries{coefficient_type(s)}(-s.coefficients)
Base.zero(s::TimeSeries) = TimeSeries{coefficient_type(s)}([zero(c) for c in s.coefficients])
Base.one(s::TimeSeries) = zero(s) + one(constant_term(s))

function time_order(a::TimeSeries, b::TimeSeries)
    max_order(a) == max_order(b) || throw(ArgumentError("Time series orders must match"))
    return max_order(a)
end
for op in (:+, :-)
    @eval function Base.$op(a::TimeSeries, b::TimeSeries)
        time_order(a, b)
        c = $op.(a.coefficients, b.coefficients)
        return TimeSeries{eltype(c)}(c)
    end
end
function Base.:+(a::TimeSeries, b::Real)
    T = promote_type(coefficient_type(a), typeof(b))
    c = T[copy(x) for x in a.coefficients]
    c[1] += b
    return TimeSeries{T}(c)
end
Base.:+(a::Real, b::TimeSeries) = b + a
Base.:-(a::TimeSeries, b::Real) = a + (-b)
Base.:-(a::Real, b::TimeSeries) = -b + a
for op in (:*, :/)
    @eval function Base.$op(a::TimeSeries, b::Real)
        c = $op.(a.coefficients, b)
        return TimeSeries{eltype(c)}(c)
    end
end
Base.:*(a::Real, b::TimeSeries) = b * a
function Base.:*(a::TimeSeries, b::TimeSeries)
    n = time_order(a, b)
    c0 = constant_term(a) * constant_term(b)
    c = Vector{typeof(c0)}(undef, n + 1)
    c[1] = c0
    @inbounds for k in 1:n
        c[k + 1] = time_convolution(a.coefficients, b.coefficients, k, c0)
    end
    return TimeSeries{eltype(c)}(c)
end
function Base.:/(a::TimeSeries, b::TimeSeries)
    n = time_order(a, b)
    b0 = constant_term(b)
    iszero(constant_term(b0)) && throw(DomainError(b0, "Division needs an invertible constant coefficient"))
    c0 = constant_term(a) / b0
    c = Vector{typeof(c0)}(undef, n + 1)
    c[1] = c0
    reciprocal = inv(b0)
    @inbounds for k in 1:n
        c[k + 1] = time_quotient_coefficient(a.coefficients, b.coefficients, c, k, reciprocal)
    end
    return TimeSeries{eltype(c)}(c)
end
Base.inv(s::TimeSeries) = one(s) / s
Base.:/(a::Real, b::TimeSeries) = a * inv(b)
function Base.:^(s::TimeSeries, n::Integer)
    # Avoid negating typemin and keep the number of products logarithmic.
    if n < 0
        inverse = inv(s)
        return inverse^(-(n + 1)) * inverse
    end
    n == 1 && return copy(s)
    result, factor = one(s), s
    while n > 0
        isodd(n) && (result = result * factor)
        n >>= 1
        n > 0 && (factor = factor * factor)
    end
    return result
end

# From a*y' = p*a'*y (power), y' = a'*y (exp), or a*y' = a'
# (log). Convolution in time leaves the coefficient algebra untouched.
function time_recurrence(a::TimeSeries, c0, p, ::Val{MODE}) where {MODE}
    n = max_order(a)
    c = Vector{typeof(c0)}(undef, n + 1)
    c[1] = c0
    a0 = constant_term(a)
    MODE !== :exp && n > 0 && iszero(constant_term(a0)) &&
        throw(DomainError(a0, "The expansion needs an invertible constant coefficient"))
    reciprocal = MODE === :exp || n == 0 ? one(c0) : inv(a0)
    @inbounds for k in 1:n
        value = time_weighted_sum(a.coefficients, c, k, p, Val(MODE))
        MODE === :log && (value += k * a.coefficients[k + 1])
        c[k + 1] = MODE === :exp ? value / k : time_normalize(value, reciprocal, k)
    end
    return TimeSeries{eltype(c)}(c)
end
Base.exp(a::TimeSeries) = time_recurrence(a, exp(constant_term(a)), 0, Val(:exp))
function Base.log(a::TimeSeries)
    a0 = constant_term(a)
    iszero(constant_term(a0)) && throw(DomainError(a0, "Logarithm needs a nonzero expansion point"))
    return time_recurrence(a, log(a0), 0, Val(:log))
end
function time_power(a::TimeSeries, p::Real)
    p isa TaylorPolynomial && return exp(p * log(a))
    isinteger(p) && return a^BigInt(p)
    a0 = constant_term(a)
    iszero(constant_term(a0)) && throw(DomainError(a0, "Fractional powers need a nonzero expansion point"))
    return time_recurrence(a, a0^p, p, Val(:power))
end
Base.:^(a::TimeSeries, p::Real) = time_power(a, p)
Base.:^(a::TimeSeries, p::Rational) = time_power(a, p)
Base.:^(a::TimeSeries, b::TimeSeries) = exp(b * log(a))
Base.:^(a::Real, b::TimeSeries) = exp(log(a) * b)
Base.:^(::Irrational{:ℯ}, b::TimeSeries) = exp(b)
Base.sqrt(a::TimeSeries) = time_recurrence(a, sqrt(constant_term(a)), 1 // 2, Val(:power))
Base.cbrt(a::TimeSeries) = time_recurrence(a, cbrt(constant_term(a)), 1 // 3, Val(:power))
function time_sincos(a::TimeSeries, hyperbolic::Bool)
    a0 = constant_term(a)
    s0, c0 = hyperbolic ? (sinh(a0), cosh(a0)) : sincos(a0)
    n = max_order(a)
    s, c = Vector{typeof(s0)}(undef, n + 1), Vector{typeof(c0)}(undef, n + 1)
    s[1], c[1] = s0, c0
    @inbounds for k in 1:n
        sv = time_weighted_sum(a.coefficients, c, k, 0, Val(:exp))
        cv = time_weighted_sum(a.coefficients, s, k, 0, Val(:exp))
        s[k + 1], c[k + 1] = sv / k, (hyperbolic ? cv : -cv) / k
    end
    return TimeSeries{eltype(s)}(s), TimeSeries{eltype(c)}(c)
end
Base.sincos(a::TimeSeries) = time_sincos(a, false)
for (fn, hyperbolic, index) in ((:sin, false, 1), (:cos, false, 2), (:sinh, true, 1), (:cosh, true, 2))
    @eval Base.$fn(a::TimeSeries) = time_sincos(a, $hyperbolic)[$index]
end
Base.tan(a::TimeSeries) = sin(a) / cos(a)
Base.tanh(a::TimeSeries) = sinh(a) / cosh(a)

function evaluate(s::TimeSeries, τ::Real)
    value = copy(last(s.coefficients))
    @inbounds for k in (length(s.coefficients) - 1):-1:1
        value = muladd(value, τ, s.coefficients[k])
    end
    return value
end
(s::TimeSeries)(τ::Real) = evaluate(s, τ)
function differentiate(s::TimeSeries)
    max_order(s) == 0 && return TimeSeries([zero(constant_term(s))])
    c = [k * s.coefficients[k + 1] for k in 1:max_order(s)]
    return TimeSeries{eltype(c)}(c)
end
function integrate(s::TimeSeries)
    c = [s.coefficients[k + 1] / (k + 1) for k in 0:max_order(s)]
    pushfirst!(c, zero(first(c)))
    return TimeSeries{eltype(c)}(c)
end

"""
    taylor_expand(f, initial, t0; order = 20, parameters = nothing)

Expand the solution of `u' = f(u, parameters, t)` about `t0`, returning one
[`TimeSeries`](@ref) per state component. `initial` is a nonempty real vector;
its entries may be multivariate Taylor polynomials. Evaluate at a **time
offset** `h` with `[s(h) for s in expansion]` to take a Taylor step.

The recurrence `u[k+1] = coefficient(f(u, parameters, t), k)/(k+1)` computes
each intermediate coefficient once. The right-hand side is evaluated first
at the expansion point to determine coefficient types, then with recording
scalars to construct its arithmetic graph.
The right-hand side must be analytic near the expansion point and return a
vector, without mutating its arguments. Use analytic arithmetic supported by
`TimeSeries`; inspecting series coefficients or branching on time or state
values inside the right-hand side is not supported.

Time order is independent of the current multivariate truncation order. This
function does not initialize or change the polynomial algebra. Numeric states
need no algebra. It constructs a local expansion, **not** an error-controlled
ODE solution; choose steps inside its domain of convergence. See the Taylor
integration tutorial for adaptive stepping and independent validation.
"""
function taylor_expand(f, initial::AbstractVector{<:Real}, t0::Real; order::Integer = 20, parameters = nothing)
    return prepare_time_expansion(f, initial, t0; order, parameters).series
end
function prepare_time_expansion(f, initial::AbstractVector{<:Real}, t0::Real; order::Integer = 20, parameters = nothing)
    order >= 1 || throw(ArgumentError("Time order must be positive"))
    isempty(initial) && throw(ArgumentError("The initial state cannot be empty"))
    isfinite(t0) && all(isfinite, initial) || throw(ArgumentError("Expansion point must be finite"))
    t0 isa TaylorPolynomial && throw(ArgumentError("The expansion epoch must be a scalar time"))
    state = [u / one(u) for u in initial]
    epoch = t0 / one(t0)
    # Evaluate the constant RHS first to determine promotion, including uncertain
    # parameters supplied with an initially numeric state.
    derivative = f(state, parameters, epoch)
    derivative isa AbstractVector && length(derivative) == length(state) ||
        throw(DimensionMismatch("The right-hand side must return one derivative per state"))
    all(x -> x isa Real, derivative) || throw(ArgumentError("Unsupported right-hand-side element"))
    all(isfinite, derivative) || throw(DomainError(derivative, "Nonfinite time coefficient"))
    T = foldl((S, x) -> promote_type(S, typeof(x)), state; init = typeof(first(state)))
    T = foldl((S, x) -> promote_type(S, typeof(x / 1)), derivative; init = T)
    return time_workspace(f, state, epoch, Int(order), parameters, T)
end

"""
    TaylorMethod(order = 20)

An adaptive Taylor ODE algorithm for the SciML `solve`/`init` interface. Load
OrdinaryDiffEqCore (or an OrdinaryDiffEq solver package) to enable this optional
extension. Supports in-place and out-of-place analytic vector ODEs, numeric or
[`TaylorPolynomial`](@ref) states, forward/backward integration, `saveat`, dense
Taylor interpolation and callbacks. Time order is independent of uncertainty
order. Choose `order ≥ 3`; the state dimension stays fixed during a solve.

The step controller uses the last two retained time coefficients, scaled by
`abstol` and `reltol`. Polynomial errors use the sum of absolute uncertainty
coefficients before applying the solver's `internalnorm`. This controls a
local truncation estimate, not a rigorous error bound. Standard SciML options
such as `dt`, `dtmax`, `adaptive=false`, and `maxiters` are supported.

```julia
using DifferentialAlgebra, OrdinaryDiffEqCore, SciMLBase
problem = ODEProblem((u, p, t) -> [u[2], -u[1]], [0.0, 1.0], (0.0, 2π))
solution = solve(problem, TaylorMethod(20); abstol = 1.0e-12, reltol = 1.0e-12)
```
"""
function TaylorMethod end
