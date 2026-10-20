mutable struct Algebra
    basis::MonomialBasis
    names::Vector{String}
    active::Bool
end
const CURRENT_ALGEBRA = Ref{Union{Nothing, Algebra}}(nothing)
const ALGEBRA_LOCK = ReentrantLock() # Serializes algebra initialization.

# Temporary contexts are strong scope bindings. Policies for the shared global
# default use weak references so an idle task cannot retain obsolete tables.
mutable struct AlgebraState
    const owner::Union{Algebra, WeakRef}
    cutoff::Int
    epsilon::Float64
    big_epsilon::Union{Nothing, BigFloat}
    const parent::Union{Nothing, AlgebraState}
end
const ALGEBRA_STATE_KEY = gensym(:DifferentialAlgebraState)
@inline function algebra_state()
    state = get(task_local_storage(), ALGEBRA_STATE_KEY, nothing)
    return state::Union{Nothing, AlgebraState}
end
@inline function state_algebra(state::AlgebraState)
    owner = state.owner
    return owner isa Algebra ? owner : owner.value::Union{Nothing, Algebra}
end
@inline function current_algebra()
    state = algebra_state()
    return state !== nothing && state.owner isa Algebra ? state.owner : CURRENT_ALGEBRA[]
end
@inline function algebra_settings(ctx::Algebra)
    state = algebra_state()
    while state !== nothing
        state_algebra(state) === ctx && return state
        state = state.parent
    end
    return nothing
end
@inline function working_order(ctx::Algebra)
    state = algebra_settings(ctx)
    return state === nothing ? ctx.basis.order : state.cutoff
end
@inline function float_tolerance(ctx::Algebra)
    state = algebra_settings(ctx)
    return state === nothing ? 0.0 : state.epsilon
end
@inline function big_tolerance(ctx::Algebra)
    state = algebra_settings(ctx)
    return state === nothing ? nothing : state.big_epsilon
end
function current_settings(ctx::Algebra)
    state = algebra_state()
    state !== nothing && state_algebra(state) === ctx && return state
    state = AlgebraState(WeakRef(ctx), ctx.basis.order, 0.0, nothing, nothing)
    task_local_storage(ALGEBRA_STATE_KEY, state)
    return state
end

# A real multivariate Taylor polynomial with coefficients of type T.
"""
    TaylorPolynomial{T}(c = zero(T))
    TaylorPolynomial(c::Real)
    TaylorPolynomial(i::Integer, c::Real)

A multivariate Taylor polynomial with coefficients of concrete real type `T`.
`TaylorPolynomial(c)` creates a constant with coefficient type `typeof(float(c))`.
`TaylorPolynomial{T}(c)` explicitly chooses storage type `T`. The two-argument constructor
creates `c` times variable `i`; index zero creates a constant.

Use [`variables`](@ref) to initialize the algebra. A polynomial is callable:
`p(point)` evaluates it numerically or composes it with polynomial coordinates.
"""
mutable struct TaylorPolynomial{T <: Real} <: Real
    coeffs::Vector{T}
    len::Int
    algebra::Algebra
    function TaylorPolynomial{T}(coeffs::Vector{T}, len::Int, algebra::Algebra) where {T <: Real}
        isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
        1 <= len <= length(coeffs) <= last(algebra.basis.ends) || throw(ArgumentError("Invalid coefficient storage length"))
        return new{T}(coeffs, len, algebra)
    end
end
TaylorPolynomial(coeffs::Vector{T}, len::Int, algebra::Algebra) where {T <: Real} = TaylorPolynomial{T}(coeffs, len, algebra)
struct TaylorError <: Exception
    message::String
end
Base.showerror(io::IO, e::TaylorError) = print(io, "DifferentialAlgebra: ", e.message)

@inline function ready()
    ctx = current_algebra()
    ctx === nothing && throw(ArgumentError("Call initialize!(order, variables) first"))
    return ctx
end
@inline function valid(a::TaylorPolynomial)
    a.algebra.active || throw(ArgumentError("Polynomial belongs to a previous initialization"))
    return a.algebra
end
@inline function compatible(a::TaylorPolynomial, b::TaylorPolynomial)
    ctx = valid(a)
    valid(b) === ctx || throw(ArgumentError("Polynomials have different contexts"))
    return ctx
end

"""
    initialize!(order, n; names = nothing, table_bytes = 32 * 1024^2)

Initialize the global default algebra with maximum degree `order` and `n` variables.
Inside a temporary construction scope, reinitialize only that task's algebra.
Existing polynomials of the replaced algebra become invalid. Both dimensions
must be positive.
Use [`variables`](@ref) to initialize and construct the variables in one call.
`table_bytes` bounds lookup tables, excluding polynomial storage and basis metadata.
Configure the algebra before launching concurrent calculations.
`names` is an optional tuple or vector of distinct symbol or string identifiers
used when displaying polynomials. The default labels are `x₁`, `x₂`, and so on.
"""
function initialize!(
        no::Integer, nv::Integer; names = nothing,
        table_bytes::Integer = 32 * 1024^2
    )
    state = algebra_state()
    if state !== nothing && state.owner isa Algebra
        old = state.owner
        ctx = make_algebra(no, nv, old; names, table_bytes)
        old.active = false
        task_local_storage(
            ALGEBRA_STATE_KEY,
            AlgebraState(ctx, ctx.basis.order, 0.0, nothing, state.parent)
        )
    else
        lock(ALGEBRA_LOCK) do
            old = CURRENT_ALGEBRA[]
            ctx = make_algebra(no, nv, old; names, table_bytes)
            old !== nothing && (old.active = false)
            CURRENT_ALGEBRA[] = ctx
            task_local_storage(
                ALGEBRA_STATE_KEY,
                AlgebraState(WeakRef(ctx), ctx.basis.order, 0.0, nothing, nothing)
            )
        end
    end
    return nothing
end
function make_algebra(
        no::Integer, nv::Integer, old;
        names = nothing, table_bytes::Integer = 32 * 1024^2
    )
    1 <= no <= 65535 && 1 <= nv <= 1024 ||
        throw(ArgumentError("Invalid order or variable count"))
    nm = binomial(big(no) + nv, nv)
    nm * nv <= 32 * 1024^2 && nm <= typemax(Int32) ||
        throw(ArgumentError("Monomial basis is too large"))
    0 <= table_bytes <= typemax(Int) ||
        throw(ArgumentError("Invalid lookup-table budget"))
    labels = variable_labels(names, Int(nv))
    reusable = old !== nothing && old.basis.order == no &&
        old.basis.variables == nv && old.basis.table_bytes == table_bytes
    basis = reusable ? old.basis :
        MonomialBasis(Int(no), Int(nv); table_bytes = Int(table_bytes))
    return Algebra(basis, labels, true)
end

function variable_labels(names, n::Int)
    names === nothing && return ["x" * join('₀' + (c - '0') for c in string(i)) for i in 1:n]
    names isa Union{Tuple, AbstractVector} || throw(ArgumentError("Names must be a tuple or vector of symbols or strings"))
    length(names) == n || throw(DimensionMismatch("Provide one name per independent variable"))
    all(name -> name isa Union{Symbol, AbstractString}, names) || throw(ArgumentError("Variable names must be symbols or strings"))
    labels = String[String(name) for name in names]
    all(Base.isidentifier, labels) || throw(ArgumentError("Variable names must be identifiers"))
    allunique(labels) || throw(ArgumentError("Variable names must be distinct"))
    return labels
end
"""
    max_order()
    max_order(p::TaylorPolynomial)

Maximum total degree of the current algebra, or of the algebra owning `p`.
For ADS maps and their patches, return the requested retained order, separately
from the actual stored `degree`. Taylor models report their fixed retained order.
"""
max_order() = ready().basis.order
# Internal initialization query for loading regressions.
isinitialized() = current_algebra() !== nothing && current_algebra().active

# Temporary maps own their compiled coefficients. Restore the caller's algebra,
# including the uninitialized state, even if the callback fails or reinitializes.
function with_algebra(f::F, order, n; kwargs...) where {F}
    ctx = make_algebra(order, n, current_algebra(); kwargs...)
    state = AlgebraState(ctx, ctx.basis.order, 0.0, nothing, algebra_state())
    return task_local_storage(ALGEBRA_STATE_KEY, state) do
        try
            return f(ctx)
        finally
            current = algebra_state()
            current !== nothing && current.owner isa Algebra &&
                (current.owner.active = false)
            ctx.active = false
        end
    end
end
"""
    nvariables()
    nvariables(p::TaylorPolynomial)
    nvariables(map::CompiledMap)

Number of independent variables in the current algebra or polynomial `p`.
For a compiled map, return the highest variable index used by its terms.
"""
nvariables() = ready().basis.variables
# Internal basis-capacity query for storage tests.
nmonomials() = last(ready().basis.ends)
"""
    truncation_order()

The current task's working order, bounded by the algebra's `max_order()`.
"""
truncation_order() = working_order(ready())
"""
    coefficient_tolerance()

The current task's absolute threshold for discarding small floating coefficients.
Returns `Float64`, or the preserved `BigFloat` value when one was explicitly set.
"""
function coefficient_tolerance()
    state = algebra_settings(ready())
    return state === nothing ? 0.0 : something(state.big_epsilon, state.epsilon)
end
"""
    set_truncation_order!(order = max_order())

Set this task's working total degree and return the previous order. The minimum working
order is one; zero selects one. Prefer `with_order` for temporary changes.
"""
function set_truncation_order!(n::Integer = max_order())
    ctx = ready()
    0 <= n <= ctx.basis.order || throw(ArgumentError("Invalid truncation order"))
    state = current_settings(ctx)
    old = state.cutoff
    state.cutoff = max(1, Int(n)) # Working orders are positive.
    return old
end
"""
    set_coefficient_tolerance!(tolerance)

Set this task's finite, nonnegative filtering threshold and return its previous
value. Exact real coefficient types retain all nonzero coefficients.
"""
function set_coefficient_tolerance!(value::Real)
    isfinite(value) && value >= 0 || throw(ArgumentError("Epsilon must be finite and nonnegative"))
    ctx = ready(); old = coefficient_tolerance()
    state = current_settings(ctx)
    state.epsilon = Float64(value)
    state.big_epsilon = value isa BigFloat ? copy(value) : nothing
    return old
end

"""
    with_order(f, order)

Run `f()` with a temporary truncation order, restoring the previous order even
if `f` throws. Use `with_order(order) do ... end`. Order and filtering settings
are task-local. Newly spawned tasks use full order and zero tolerance unless
they configure their own settings; temporary construction scopes are not inherited.
Arithmetic results, including scalar addition, unary signs and integer powers,
discard terms above the working order. Copying, coefficient conversion and
inspection preserve the stored polynomial; changing the order alone does not
alter existing storage. Discarding terms is not a certified error estimate.
"""
function with_order(f::F, order::Integer) where {F}
    state = current_settings(ready())
    previous = set_truncation_order!(order)
    try
        return f()
    finally
        state.cutoff = previous
    end
end

max_order(p::TaylorPolynomial) = valid(p).basis.order
nvariables(p::TaylorPolynomial) = valid(p).basis.variables

# Narrow scalar hooks for optional coefficient types. Ordinary real kernels
# inline to their original operations; extensions must preserve scalar metadata.
@inline coefficient_iszero(c::Real) = iszero(c)
@inline coefficient_isone(c::Real) = isone(c)
@inline coefficient_isinteger(c::Real) = isinteger(c)
@inline coefficient_isequal(a::Real, b::Real) = isequal(a, b)
@inline coefficient_convert(::Type{T}, c::Real) where {T <: Real} = convert(T, c)
@inline coefficient_operand(::Type{T}, c::Real) where {T <: Real} = c
@inline degree_factor(c::Real, n::Integer) = n
@inline coefficient_muladd(a::Real, b::Real, c::Real) = muladd(a, b, c)
@inline coefficient_power(c::Real, n::Integer) = c^n
@inline invertible_constant(c::Real) = !coefficient_iszero(c)

@inline keep(c::Real, epsilon::Float64) = !coefficient_iszero(c)
@inline keep(c::AbstractFloat, epsilon::Float64) = !(abs(c) <= epsilon)
@inline keep(c::Real, epsilon::Float64, precise) = keep(c, epsilon)
@inline keep(c::AbstractFloat, epsilon::Float64, precise) = precise === nothing ?
    keep(c, epsilon) : !(abs(c) <= precise)
@inline keep(c::BigFloat, epsilon::Float64, precise) =
    !(abs(c) <= something(precise, epsilon))
@inline function keep(c::Real, ctx::Algebra)
    state = algebra_settings(ctx)
    return state === nothing ? keep(c, 0.0) :
        keep(c, state.epsilon, state.big_epsilon)
end
@inline function finish!(a::TaylorPolynomial, n::Int = length(a.coeffs); filter::Bool = true)
    if filter
        state = algebra_settings(a.algebra)
        epsilon = state === nothing ? 0.0 : state.epsilon
        precise = state === nothing ? nothing : state.big_epsilon
        if epsilon != 0 || precise !== nothing
            @inbounds for i in 1:n
                keep(a.coeffs[i], epsilon, precise) ||
                    (a.coeffs[i] = zero(eltype(a.coeffs)))
            end
        end
    end
    @inbounds while n > 1 && coefficient_iszero(a.coeffs[n])
        n -= 1
    end
    a.len = n
    return a
end
@inline allocate(ctx::Algebra, ::Type{T}, n::Int = last(ctx.basis.ends)) where {T <: Real} = TaylorPolynomial{T}(zeros(T, n), 1, ctx)
# Only for kernels that write every active coefficient before returning. As with
# a reused buffer, capacity beyond len is unspecified and must never be read.
@inline allocate_undef(ctx::Algebra, ::Type{T}, n::Int) where {T <: Real} = TaylorPolynomial{T}(Vector{T}(undef, n), 1, ctx)
@inline function reserve!(a::TaylorPolynomial, n::Int)
    old = length(a.coeffs)
    if n > old
        resize!(a.coeffs, n)
        @inbounds for i in (old + 1):n
            a.coeffs[i] = zero(eltype(a.coeffs))
        end
    end
    return a
end
"""
    coefficient_type(p)
    coefficient_type(TaylorPolynomial{T})

Return the scalar coefficient type `T`.
"""
coefficient_type(::TaylorPolynomial{T}) where {T} = T
coefficient_type(::Type{TaylorPolynomial{T}}) where {T} = T

TaylorPolynomial() = TaylorPolynomial(0.0)
TaylorPolynomial(x::Real) = TaylorPolynomial{typeof(float(x))}(x)
TaylorPolynomial(a::TaylorPolynomial) = copy(a)
function TaylorPolynomial{T}() where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    return TaylorPolynomial{T}(zero(T))
end
function TaylorPolynomial{T}(x::Real) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    return TaylorPolynomial{T}(ready(), x)
end
# Explicit owning contexts keep object-based operations independent of factories
# bound to the current task's temporary algebra.
function TaylorPolynomial{T}(ctx::Algebra, x::Real) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    c = coefficient_convert(T, x)
    return TaylorPolynomial{T}([keep(c, ctx) ? c : zero(T)], 1, ctx)
end
function TaylorPolynomial{T}(i::Integer, c::Real) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    return TaylorPolynomial{T}(ready(), i, c)
end
function TaylorPolynomial{T}(ctx::Algebra, i::Integer, c::Real) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    0 <= i <= ctx.basis.variables || throw(ArgumentError("Variable index out of bounds"))
    a = allocate(ctx, T, Int(i) + 1)
    a.coeffs[i + 1] = coefficient_convert(T, c)
    return finish!(a)
end
TaylorPolynomial(i::Integer, c::Real) = TaylorPolynomial{typeof(float(c))}(i, c)
"""
    variable(i, T = Float64)

Return independent variable `i` in the current algebra with coefficient type `T`.
Indices start at one. This does not reinitialize or invalidate other polynomials.
"""
function variable(i::Integer, ::Type{T} = Float64) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    i >= 1 || throw(ArgumentError("Variable indices start at one"))
    return TaylorPolynomial{T}(i, one(T))
end
"""
    constant_term(p)

Return the constant coefficient of a polynomial, or the constant parts of an array.
For a real scalar, return the scalar itself.
"""
@inline constant_term(a::TaylorPolynomial) = (valid(a); @inbounds a.coeffs[1])
constant_term(a::Real) = a
constant_term(a::AbstractArray{<:TaylorPolynomial}) = constant_term.(a)
constant_term(a::AbstractArray{<:Real}) = a
"""
    linear_part(p)
    linear_part(map)

Extract first-degree coefficients as a vector, or one row per map component.
This equals the gradient or Jacobian evaluated at the expansion origin.
"""
function linear_part(a::TaylorPolynomial{T}) where {T}
    ctx = valid(a)
    out = zeros(T, ctx.basis.variables)
    @inbounds for i in 1:min(length(out), a.len - 1)
        out[i] = a.coeffs[i + 1]
    end
    return out
end
function linear_part(v::AbstractVector{<:TaylorPolynomial})
    T = isempty(v) ? Float64 : mapreduce(coefficient_type, promote_type, v)
    n = isempty(v) ? nvariables() : nvariables(first(v))
    result = zeros(T, length(v), n)
    for (i, p) in enumerate(v)
        result[i, :] = linear_part(p)
    end
    return result
end

"""
    variables(n; order, names = nothing, table_bytes = 32 * 1024^2)
    variables(T, n; order, names = nothing, table_bytes = 32 * 1024^2)
    variables(names; order, table_bytes = 32 * 1024^2)
    variables(T, names; order, table_bytes = 32 * 1024^2)

Initialize an algebra of total degree `order` and return its `n` independent variables.
`T` is a concrete real coefficient type and defaults to `Float64`. The result is
an ordinary `Vector{TaylorPolynomial{T}}`. Reinitialization invalidates existing polynomials;
use [`variable`](@ref) to retrieve a variable in the current algebra.

Supply a tuple or vector of distinct symbol or string identifiers to name the
variables in polynomial displays. The number of variables can be inferred from
the names, or supplied explicitly with the `names` keyword. Names are copied at
initialization and do not affect arithmetic or indexing. Defaults are `x₁`, `x₂`, …;
custom names are displayed verbatim, and powers use Unicode superscripts.

# Examples
```julia
x, y = variables((:x, :y); order = 6)
p = sin(x) * exp(y)
p([0.1, 0.2])
```
"""
function variables(::Type{T}, n::Integer; order::Integer, names = nothing, table_bytes::Integer = 32 * 1024^2) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    initialize!(order, n; names, table_bytes)
    return [variable(i, T) for i in 1:n]
end
variables(n::Integer; kwargs...) = variables(Float64, n; kwargs...)
function variables(
        ::Type{T}, names::Union{Tuple, AbstractVector}; kwargs...
    ) where {T <: Real}
    haskey(kwargs, :names) &&
        throw(ArgumentError("Supply variable names positionally or by keyword, not both"))
    return variables(T, length(names); names, kwargs...)
end
variables(names::Union{Tuple, AbstractVector}; kwargs...) = variables(Float64, names; kwargs...)
