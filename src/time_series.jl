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
Division requires an invertible denominator center. Noninteger powers and
roots require a regular center, except that a known identically zero series
has zero positive roots. A zero primal with nonzero dual partials is singular.
"""
struct TimeSeries{T <: Real} <: Number
    coefficients::Vector{T}
    function TimeSeries{T}(coefficients::Vector{T}) where {T <: Real}
        isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
        isempty(coefficients) && throw(ArgumentError("A time series needs a constant coefficient"))
        return new{T}(coefficients)
    end
end
TimeSeries(c::AbstractVector{T}) where {T <: Real} = TimeSeries{T}(T[copy(x) for x in c])
time_finish_coefficient(x::Real) = x
time_finish_coefficient(x::TaylorPolynomial) =
    finish!(x, min(x.len, x.algebra.basis.ends[working_order(x.algebra) + 1]))
function time_finish(series::TimeSeries, filter::Bool = true)
    filter && foreach(time_finish_coefficient, series.coefficients)
    return series
end
time_result(c::Vector; filter::Bool = true) = time_finish(TimeSeries{eltype(c)}(c), filter)
time_identity(s::TimeSeries) = time_add(zero(s), time_identity(constant_term(s)); filter = false)
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
Base.:+(s::TimeSeries) = time_finish(copy(s))
Base.:-(s::TimeSeries) = time_result(time_coefficient.(-, s.coefficients))
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
function time_add(a::TimeSeries, b::Real; filter::Bool = true)
    T = promote_type(coefficient_type(a), typeof(b))
    c = T[copy(x) for x in a.coefficients]
    c[1] = time_coefficient(+, c[1], b)
    return time_finish(TimeSeries{T}(c), filter)
end
Base.:+(a::TimeSeries, b::Real) = time_add(a, b)
Base.:+(a::Real, b::TimeSeries) = b + a
Base.:-(a::TimeSeries, b::Real) = a + (-b)
Base.:-(a::Real, b::TimeSeries) = -b + a
function time_multiply(a::TimeSeries, b::Real; filter::Bool = true)
    return time_result([time_coefficient(*, c, b) for c in a.coefficients]; filter)
end
Base.:*(a::TimeSeries, b::Real) = time_multiply(a, b)
function time_divide(a::TimeSeries, b::Real; filter::Bool = true)
    invertible_constant(constant_term(b)) ||
        throw(TaylorError("Division needs an invertible constant coefficient"))
    return time_result([time_coefficient(/, c, b) for c in a.coefficients]; filter)
end
Base.:/(a::TimeSeries, b::Real) = time_divide(a, b)
Base.:*(a::Real, b::TimeSeries) = b * a
function time_multiply(a::TimeSeries, b::TimeSeries; filter::Bool = true)
    n = time_order(a, b)
    c0 = time_coefficient(*, constant_term(a), constant_term(b))
    c = Vector{typeof(c0)}(undef, n + 1)
    c[1] = c0
    @inbounds for k in 1:n
        c[k + 1] = time_convolution(a.coefficients, b.coefficients, k, c0)
    end
    return time_result(c; filter)
end
Base.:*(a::TimeSeries, b::TimeSeries) = time_multiply(a, b)
function time_divide(a::TimeSeries, b::TimeSeries; filter::Bool = true)
    n = time_order(a, b)
    b0 = constant_term(b)
    !invertible_constant(constant_term(b0)) && throw(TaylorError("Division needs an invertible constant coefficient"))
    c0 = time_coefficient(/, constant_term(a), b0)
    c = Vector{typeof(c0)}(undef, n + 1)
    c[1] = c0
    @inbounds for k in 1:n
        c[k + 1] = time_quotient_coefficient(a.coefficients, b.coefficients, c, k, b0)
    end
    return time_result(c; filter)
end
Base.:/(a::TimeSeries, b::TimeSeries) = time_divide(a, b)
Base.inv(s::TimeSeries) = time_divide(time_identity(s), s)
Base.:/(a::Real, b::TimeSeries) = time_divide(time_add(zero(b), a; filter = false), b)
function time_power(s::TimeSeries, n::Integer; filter::Bool = true)
    # Avoid negating typemin and keep the number of products logarithmic.
    if n < 0
        inverse = time_divide(time_identity(s), s; filter = false)
        return time_multiply(time_power(inverse, -(n + 1); filter = false), inverse; filter)
    end
    n == 1 && return time_finish(copy(s), filter)
    result, factor = time_identity(s), s
    while n > 0
        isodd(n) && (result = time_multiply(result, factor; filter = false))
        n >>= 1
        n > 0 && (factor = time_multiply(factor, factor; filter = false))
    end
    return time_finish(result, filter)
end
Base.:^(s::TimeSeries, n::Integer) = time_power(s, n)
# Base's literal negative-power fallback would finish inv(s) too early.
Base.literal_pow(::typeof(^), s::TimeSeries, ::Val{N}) where {N} = s^N

# From a*y' = p*a'*y (power), y' = a'*y (exp), or a*y' = a'
# (log). Convolution in time leaves the coefficient algebra untouched.
function time_recurrence(a::TimeSeries, c0, p, ::Val{MODE}; filter::Bool = true) where {MODE}
    MODE === :power && iszero(a) && p > 0 && return zero(a) + zero(c0)
    n = max_order(a)
    c = Vector{typeof(c0)}(undef, n + 1)
    a0 = constant_term(a)
    c[1] = c0
    # A time-constant dual or polynomial can still carry coefficient perturbations.
    MODE !== :exp &&
        !invertible_constant(constant_term(a0)) &&
        throw(TaylorError("The expansion needs an invertible constant coefficient"))
    exponent = MODE === :power ? time_exponent(p, c0) : p
    @inbounds for k in 1:n
        c[k + 1] = time_unary_coefficient(
            a.coefficients, c, k, exponent, a0, Val(MODE)
        )
    end
    return time_result(c; filter)
end
time_exponential(a::TimeSeries; filter::Bool = true) =
    time_recurrence(a, time_coefficient(exp, constant_term(a)), 0, Val(:exp); filter)
Base.exp(a::TimeSeries) = time_exponential(a)
function time_logarithm(a::TimeSeries; filter::Bool = true)
    a0 = constant_term(a)
    !invertible_constant(constant_term(a0)) && throw(TaylorError("Logarithm needs an invertible expansion point"))
    return time_recurrence(a, domain_call(x -> time_coefficient(log, x), a0), 0, Val(:log); filter)
end
Base.log(a::TimeSeries) = time_logarithm(a)
function time_power(a::TimeSeries, p::Real; filter::Bool = true)
    p isa TaylorPolynomial && return time_exponential(time_multiply(time_logarithm(a; filter = false), p; filter = false); filter)
    p isa Union{AbstractFloat, Rational} && isinteger(p) && return time_power(a, BigInt(p); filter)
    a0 = constant_term(a)
    exponent = time_exponent(p, a0)
    return time_recurrence(a, domain_call(x -> time_coefficient(^, x, exponent), a0), exponent, Val(:power); filter)
end
Base.:^(a::TimeSeries, p::Real) = time_power(a, p)
Base.:^(a::TimeSeries, p::Rational) = time_power(a, p)
function Base.:^(a::TimeSeries, b::TimeSeries)
    time_order(a, b)
    all(iszero, @view(b.coefficients[2:end])) && return time_power(a, constant_term(b))
    return time_exponential(time_multiply(b, time_logarithm(a; filter = false); filter = false))
end
function Base.:^(a::Real, b::TimeSeries)
    R = promote_type(typeof(a), typeof(constant_term(constant_term(b))))
    base = coefficient_convert(R, a)
    coefficient_iszero(base) &&
        return zero(b) + zero_base_power(base, constant_term(b), iszero(b))
    all(iszero, @view(b.coefficients[2:end])) &&
        return time_add(zero(b), domain_call(x -> time_coefficient(^, float(base), x), constant_term(b)))
    return time_exponential(time_multiply(b, domain_call(log, base); filter = false))
end
Base.:^(::Irrational{:ℯ}, b::TimeSeries) = exp(b)
Base.sqrt(a::TimeSeries) = time_recurrence(
    a, domain_call(x -> time_coefficient(sqrt, x), constant_term(a)), 1 // 2, Val(:power)
)
Base.cbrt(a::TimeSeries) = time_recurrence(
    a, domain_call(x -> time_coefficient(cbrt, x), constant_term(a)), 1 // 3, Val(:power)
)
function time_sincos(a::TimeSeries, hyperbolic::Bool; filter::Bool = true)
    a0 = constant_term(a)
    s0, c0 = hyperbolic ? (time_coefficient(sinh, a0), time_coefficient(cosh, a0)) : time_coefficient(sincos, a0)
    n = max_order(a)
    s, c = Vector{typeof(s0)}(undef, n + 1), Vector{typeof(c0)}(undef, n + 1)
    s[1], c[1] = s0, c0
    @inbounds for k in 1:n
        sv = time_weighted_sum(a.coefficients, c, k, 0, Val(:exp))
        cv = time_weighted_sum(a.coefficients, s, k, 0, Val(:exp))
        s[k + 1], c[k + 1] = sv, hyperbolic ? cv : time_coefficient(-, cv)
    end
    return time_result(s; filter), time_result(c; filter)
end
Base.sincos(a::TimeSeries) = time_sincos(a, false)
for (fn, hyperbolic, index) in ((:sin, false, 1), (:cos, false, 2), (:sinh, true, 1), (:cosh, true, 2))
    @eval Base.$fn(a::TimeSeries) = time_sincos(a, $hyperbolic)[$index]
end
function Base.tan(a::TimeSeries)
    s, c = time_sincos(a, false; filter = false)
    return time_divide(s, c)
end
function Base.tanh(a::TimeSeries)
    s, c = time_sincos(a, true; filter = false)
    return time_divide(s, c)
end

function evaluate(s::TimeSeries, τ::Real)
    value = copy(last(s.coefficients))
    @inbounds for k in (length(s.coefficients) - 1):-1:1
        value = time_muladd(value, τ, s.coefficients[k])
    end
    return time_finish_coefficient(value)
end
(s::TimeSeries)(τ::Real) = evaluate(s, τ)
function differentiate(s::TimeSeries)
    max_order(s) == 0 && return TimeSeries([zero(constant_term(s))])
    c = [time_coefficient(*, degree_factor(s.coefficients[k + 1], k), s.coefficients[k + 1]) for k in 1:max_order(s)]
    return time_result(c)
end
function integrate(s::TimeSeries)
    c = [time_coefficient(/, s.coefficients[k + 1], degree_factor(s.coefficients[k + 1], k + 1)) for k in 0:max_order(s)]
    pushfirst!(c, zero(first(c)))
    return time_result(c)
end

"""
    taylor_expand(f, initial, t0; order = 20, parameters = nothing)

Expand the solution of `u' = f(u, parameters, t)` about `t0`, returning one
[`TimeSeries`](@ref) per state component. `initial` is a nonempty real vector;
its entries may be multivariate Taylor polynomials. Evaluate at a **time
offset** `h` with `[s(h) for s in expansion]` to take a Taylor step.

The recurrence `u[k+1] = coefficient(f(u, parameters, t), k)/(k+1)` computes
each intermediate coefficient once. The right-hand side is evaluated twice
with recording scalars: first to determine coefficient types at the expansion
point, then to construct its concrete arithmetic graph.
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
    Base.require_one_based_indexing(initial)
    1 <= order < typemax(Int) || throw(ArgumentError("Time order must be a positive machine integer"))
    isempty(initial) && throw(ArgumentError("The initial state cannot be empty"))
    isfinite(t0) && all(isfinite, initial) || throw(ArgumentError("Expansion point must be finite"))
    t0 isa TaylorPolynomial && throw(ArgumentError("The expansion epoch must be a scalar time"))
    state = [time_coefficient(/, u, time_identity(u)) for u in initial]
    epoch = t0 / one(t0)
    # Evaluate the constant RHS first to determine promotion, including uncertain
    # parameters supplied with an initially numeric state.
    derivative = time_initial_derivative(f, state, epoch, parameters)
    derivative isa AbstractVector && length(derivative) == length(state) ||
        throw(DimensionMismatch("The right-hand side must return one derivative per state"))
    Base.require_one_based_indexing(derivative)
    all(x -> x isa Real, derivative) || throw(ArgumentError("Unsupported right-hand-side element"))
    all(isfinite, derivative) || throw(TaylorError("Nonfinite time coefficient"))
    T = foldl((S, x) -> promote_type(S, typeof(x)), state; init = typeof(first(state)))
    T = foldl((S, x) -> promote_type(S, typeof(time_coefficient(/, x, 1))), derivative; init = T)
    return time_workspace(f, state, epoch, Int(order), parameters, T)
end

"""
    TaylorMethod(order = 20)

An adaptive Taylor ODE algorithm for the SciML `solve`/`init` interface. Load
OrdinaryDiffEqCore (or an OrdinaryDiffEq solver package) to enable this optional
extension. Supports in-place and out-of-place analytic vector ODEs, numeric or
[`TaylorPolynomial`](@ref) states, forward/backward integration, `saveat`, dense
Taylor interpolation and callbacks. Time order is independent of uncertainty
order. Choose `order ≥ 3`; use a dense `Vector` whose dimension stays fixed during a solve.
Use polynomial initial states for uncertain polynomial parameters. Interval
states and parameters are rejected at initialization; this algorithm does not
certify time-integration errors.

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
TaylorMethod() = TaylorMethod(20)
TaylorMethod(order::Integer, args...; kwargs...) =
    throw(ArgumentError("Load OrdinaryDiffEqCore to use TaylorMethod"))

# Optional scalar extensions reject unsupported controller scalars without
# coupling the solver extension to those dependencies.
taylor_solver_coefficient(::Type{<:Real}) = true
taylor_solver_coefficient(::Type{TaylorPolynomial{T}}) where {T} =
    taylor_solver_coefficient(T)
