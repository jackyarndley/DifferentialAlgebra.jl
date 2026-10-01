mutable struct Algebra
    basis::MonomialBasis
    names::Vector{String}
    cutoff::Int
    epsilon::Float64
    big_epsilon::Union{Nothing, BigFloat}
    active::Bool
    stack::Vector{Int}
end
const CURRENT_ALGEBRA = Ref{Union{Nothing, Algebra}}(nothing)
const ALGEBRA_LOCK = ReentrantLock() # Serializes algebra initialization.

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
end
struct TaylorError <: Exception
    message::String
end
Base.showerror(io::IO, e::TaylorError) = print(io, "DifferentialAlgebra: ", e.message)

@inline function ready()
    ctx = CURRENT_ALGEBRA[]
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

Initialize the global algebra with maximum total degree `order` and `n` variables.
Existing polynomials become invalid. Both dimensions must be positive.
Use [`variables`](@ref) to initialize and construct the variables in one call.
`table_bytes` bounds lookup tables, excluding polynomial storage and basis metadata.
Configure the algebra before launching concurrent calculations.
`names` is an optional tuple or vector of distinct symbol or string identifiers
used when displaying polynomials. The default labels are `x₁`, `x₂`, and so on.
"""
function initialize!(no::Integer, nv::Integer; names = nothing, table_bytes::Integer = 32 * 1024^2)
    1 <= no <= 65535 && 1 <= nv <= 1024 || throw(ArgumentError("Invalid order or variable count"))
    nm = binomial(big(no) + nv, nv)
    nm * nv <= 32 * 1024^2 && nm <= typemax(Int32) || throw(ArgumentError("Monomial basis is too large"))
    0 <= table_bytes <= typemax(Int) || throw(ArgumentError("Invalid lookup-table budget"))
    labels = variable_labels(names, Int(nv))
    lock(ALGEBRA_LOCK) do
        old = CURRENT_ALGEBRA[]
        basis = old !== nothing && old.basis.order == no && old.basis.variables == nv && old.basis.table_bytes == table_bytes ? old.basis : MonomialBasis(Int(no), Int(nv); table_bytes = Int(table_bytes))
        old !== nothing && (old.active = false)
        CURRENT_ALGEBRA[] = Algebra(basis, labels, Int(no), 0.0, nothing, true, Int[])
    end
    return nothing
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
"""
max_order() = ready().basis.order
isinitialized() = CURRENT_ALGEBRA[] !== nothing && CURRENT_ALGEBRA[].active

# Temporary maps own their compiled coefficients. Restore the caller's algebra,
# including the uninitialized state, even if the callback fails or reinitializes.
function with_algebra(f, order, n; kwargs...)
    return lock(ALGEBRA_LOCK) do
        previous = CURRENT_ALGEBRA[]
        active = previous !== nothing && previous.active
        temporary = nothing
        try
            initialize!(order, n; kwargs...)
            temporary = ready()
            return f(temporary)
        finally
            current = CURRENT_ALGEBRA[]
            current !== nothing && current !== previous && (current.active = false)
            temporary !== nothing && (temporary.active = false)
            CURRENT_ALGEBRA[] = previous
            previous !== nothing && (previous.active = active)
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
nmonomials() = last(ready().basis.ends)
"""
    truncation_order()

Current working order, bounded by the algebra's `max_order()`.
"""
truncation_order() = ready().cutoff
"""
    coefficient_tolerance()

Absolute threshold used to discard small floating-point coefficients.
"""
coefficient_tolerance() = something(ready().big_epsilon, ready().epsilon)
"""
    set_truncation_order!(order = max_order())

Set the working total degree and return the previous order. The minimum working
order is one; zero selects one. Prefer `with_order` for temporary changes.
"""
function set_truncation_order!(n::Integer = max_order())
    ctx = ready()
    0 <= n <= ctx.basis.order || throw(ArgumentError("Invalid truncation order"))
    old = ctx.cutoff
    ctx.cutoff = max(1, Int(n)) # Working orders are positive.
    return old
end
"""
    set_coefficient_tolerance!(tolerance)

Set a finite, nonnegative coefficient filtering threshold and return its previous
value. Exact real coefficient types retain all nonzero coefficients.
"""
function set_coefficient_tolerance!(value::Real)
    isfinite(value) && value >= 0 || throw(ArgumentError("Epsilon must be finite and nonnegative"))
    ctx = ready(); old = coefficient_tolerance()
    ctx.epsilon = Float64(value)
    ctx.big_epsilon = value isa BigFloat ? copy(value) : nothing
    return old
end
push_order!(n::Integer = max_order()) = (push!(ready().stack, set_truncation_order!(n)); nothing)
function pop_order!()
    isempty(ready().stack) && throw(ArgumentError("Truncation-order stack is empty"))
    set_truncation_order!(pop!(ready().stack))
    return nothing
end

"""
    with_order(f, order)

Run `f()` with a temporary truncation order, restoring the previous order even
if `f` throws. Use `with_order(order) do ... end`. Like other algebra settings,
this is global configuration and must not change during concurrent calculations.
"""
function with_order(f::F, order::Integer) where {F}
    algebra = ready()
    previous = set_truncation_order!(order)
    try
        return f()
    finally
        algebra.cutoff = previous
    end
end

max_order(p::TaylorPolynomial) = valid(p).basis.order
nvariables(p::TaylorPolynomial) = valid(p).basis.variables

@inline keep(c::Real, epsilon::Float64) = !iszero(c)
@inline keep(c::AbstractFloat, epsilon::Float64) = !(abs(c) <= epsilon)
@inline keep(c::Real, ctx::Algebra) = keep(c, ctx.epsilon)
@inline keep(c::AbstractFloat, ctx::Algebra) = ctx.big_epsilon === nothing ? keep(c, ctx.epsilon) : !(abs(c) <= ctx.big_epsilon)
@inline keep(c::BigFloat, ctx::Algebra) = !(abs(c) <= something(ctx.big_epsilon, ctx.epsilon))
@inline function finish!(a::TaylorPolynomial, n::Int = length(a.coeffs); filter::Bool = true)
    if filter && (a.algebra.epsilon != 0 || a.algebra.big_epsilon !== nothing)
        @inbounds for i in 1:n
            keep(a.coeffs[i], a.algebra) || (a.coeffs[i] = zero(eltype(a.coeffs)))
        end
    end
    @inbounds while n > 1 && iszero(a.coeffs[n])
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
function TaylorPolynomial{T}(x::Real = zero(T)) where {T <: Real}
    ctx = ready()
    c = convert(T, x)
    return TaylorPolynomial{T}([keep(c, ctx) ? c : zero(T)], 1, ctx)
end
function TaylorPolynomial{T}(i::Integer, c::Real) where {T <: Real}
    ctx = ready()
    0 <= i <= ctx.basis.variables || throw(ArgumentError("Variable index out of bounds"))
    a = allocate(ctx, T, Int(i) + 1)
    a.coeffs[i + 1] = convert(T, c)
    return finish!(a)
end
TaylorPolynomial(i::Integer, c::Real) = TaylorPolynomial{typeof(float(c))}(i, c)
"""
    variable(i, T = Float64)

Return independent variable `i` in the current algebra with coefficient type `T`.
Indices start at one. This does not reinitialize or invalidate other polynomials.
"""
function variable(i::Integer, ::Type{T} = Float64) where {T <: Real}
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
    result = zeros(T, length(v), nvariables())
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
variables(::Type{T}, names::Union{Tuple, AbstractVector}; kwargs...) where {T <: Real} = variables(T, length(names); names, kwargs...)
variables(names::Union{Tuple, AbstractVector}; kwargs...) = variables(Float64, names; kwargs...)
