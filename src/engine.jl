"Degree-ordered monomials and a bounded, packed multiplication lookup table."
struct Basis
    order::Int
    variables::Int
    table_bytes::Int
    exponents::Matrix{Int}
    degrees::Vector{Int}
    ends::Vector{Int}
    counts::Matrix{Int}
    offsets::Vector{Int}
    products::Vector{Int32}
    splits::Vector{Int}
    left::Vector{Int32}
    right::Vector{Int32}
    parents::Vector{Int}
    factors::Vector{Int}
    traversal::Vector{Int}
end

# Rank a reverse-lexicographic monomial without allocating a tuple or dictionary.
@inline function rank(b::Basis, exponents, degree::Int)
    index = degree == 0 ? 1 : b.ends[degree] + 1
    remaining = degree
    @inbounds for v in 1:(b.variables - 1)
        remaining -= exponents[v]
        remaining > 0 && (index += b.counts[remaining, b.variables - v + 1])
    end
    return index
end
@inline function product_rank(b::Basis, i::Int, j::Int)
    degree = b.degrees[i] + b.degrees[j]
    index = degree == 0 ? 1 : b.ends[degree] + 1
    @inbounds for v in 1:(b.variables - 1)
        degree -= b.exponents[v, i] + b.exponents[v, j]
        degree > 0 && (index += b.counts[degree, b.variables - v + 1])
    end
    return index
end

function foreach_monomial(f, no::Int, nv::Int)
    current = zeros(Int, nv)
    function visit(v, remaining, degree)
        return if v == nv
            current[v] = remaining
            f(current, degree)
        else
            for value in remaining:-1:0
                current[v] = value
                visit(v + 1, remaining - value, degree)
            end
        end
    end
    for degree in 0:no
        visit(1, degree, degree)
    end
    return nothing
end

function Basis(no::Int, nv::Int; table_bytes::Int = 32 * 1024^2)
    ends = [Int(binomial(big(nv + d), nv)) for d in 0:no]
    nm = last(ends)
    counts = [Int(binomial(big(d + v), v)) for d in 0:no, v in 0:nv]
    exponents = Matrix{Int}(undef, nv, nm)
    degrees = Vector{Int}(undef, nm)
    column = 0
    foreach_monomial(no, nv) do powers, degree
        column += 1
        exponents[:, column] = powers
        degrees[column] = degree
    end
    current = zeros(Int, nv)
    offsets, products = zeros(Int, nm), Int32[]
    parents, factors = ones(Int, nm), zeros(Int, nm)
    traversal = Int[]
    splits, left, right = zeros(Int, nm + 1), Int32[], Int32[]
    b = Basis(no, nv, table_bytes, exponents, degrees, ends, counts, offsets, products, splits, left, right, parents, factors, traversal)
    pairs = sum(i -> ends[no - degrees[i] + 1], 1:nm)
    if pairs <= table_bytes ÷ (3sizeof(Int32))
        sizehint!(products, pairs)
        for i in 1:nm
            offsets[i] = length(products)
            for j in 1:ends[no - degrees[i] + 1]
                push!(products, Int32(product_rank(b, i, j)))
            end
        end
        # Group nonconstant-left products by output monomial for recurrences.
        for i in 2:nm, j in 1:ends[no - degrees[i] + 1]
            k = Int(products[offsets[i] + j])
            splits[k + 1] += 1
        end
        cumsum!(splits, splits)
        resize!(left, last(splits)); resize!(right, last(splits))
        positions = copy(splits)
        for i in 2:nm, j in 1:ends[no - degrees[i] + 1]
            k = Int(products[offsets[i] + j])
            positions[k] += 1
            left[positions[k]], right[positions[k]] = i, j
        end
    end
    children = [Int[] for _ in 1:nm]
    for i in 2:nm
        current .= @view exponents[:, i]
        v = findlast(!iszero, current)::Int
        current[v] -= 1
        parents[i], factors[i] = rank(b, current, degrees[i] - 1), v
        push!(children[parents[i]], i)
    end
    stack = [1]
    while !isempty(stack)
        i = pop!(stack)
        push!(traversal, i)
        append!(stack, Iterators.reverse(children[i]))
    end
    return b
end

mutable struct Context
    basis::Basis
    cutoff::Int
    epsilon::Float64
    big_epsilon::Union{Nothing, BigFloat}
    active::Bool
    stack::Vector{Int}
end
const engine = Ref{Context}()
const engine_lock = ReentrantLock() # Serializes algebra initialization.

# A real multivariate Taylor polynomial with coefficients of type T.
mutable struct DA{T <: Real} <: Real
    coeffs::Vector{T}
    len::Int
    context::Context
end
struct DAError <: Exception
    message::String
end
Base.showerror(io::IO, e::DAError) = print(io, "DifferentialAlgebra: ", e.message)

@inline function ready()
    isassigned(engine) || throw(ArgumentError("Call init(order, variables) first"))
    return engine[]
end
@inline function valid(a::DA)
    a.context.active || throw(ArgumentError("Polynomial belongs to a previous initialization"))
    return a.context
end
@inline function compatible(a::DA, b::DA)
    ctx = valid(a)
    valid(b) === ctx || throw(ArgumentError("Polynomials have different contexts"))
    return ctx
end

function init(no::Integer, nv::Integer; table_bytes::Integer = 32 * 1024^2)
    1 <= no <= 65535 && 1 <= nv <= 1024 || throw(ArgumentError("Invalid order or variable count"))
    nm = binomial(big(no) + nv, nv)
    nm * nv <= 32 * 1024^2 && nm <= typemax(Int32) || throw(ArgumentError("Monomial basis is too large"))
    0 <= table_bytes <= typemax(Int) || throw(ArgumentError("Invalid lookup-table budget"))
    lock(engine_lock) do
        old = isassigned(engine) ? engine[] : nothing
        basis = old !== nothing && old.basis.order == no && old.basis.variables == nv && old.basis.table_bytes == table_bytes ? old.basis : Basis(Int(no), Int(nv); table_bytes = Int(table_bytes))
        old !== nothing && (old.active = false)
        engine[] = Context(basis, Int(no), 0.0, nothing, true, Int[])
    end
    return nothing
end
getMaxOrder() = ready().basis.order
isInitialized() = isassigned(engine) && engine[].active
getMaxVariables() = ready().basis.variables
getMaxMonomials() = last(ready().basis.ends)
getTO() = ready().cutoff
getEps() = something(ready().big_epsilon, ready().epsilon)
getEpsMac(::Type{T} = Float64) where {T <: AbstractFloat} = eps(T)
function setTO(n::Integer = getMaxOrder())
    ctx = ready()
    0 <= n <= ctx.basis.order || throw(ArgumentError("Invalid truncation order"))
    old = ctx.cutoff
    ctx.cutoff = max(1, Int(n)) # preserve DACE's legacy minimum order
    return old
end
function setEps(value::Real)
    isfinite(value) && value >= 0 || throw(ArgumentError("Epsilon must be finite and nonnegative"))
    ctx = ready(); old = getEps()
    ctx.epsilon = Float64(value)
    ctx.big_epsilon = value isa BigFloat ? copy(value) : nothing
    return old
end
pushTO(n::Integer = getMaxOrder()) = (push!(ready().stack, setTO(n)); nothing)
function popTO()
    isempty(ready().stack) && throw(ArgumentError("Truncation-order stack is empty"))
    setTO(pop!(ready().stack))
    return nothing
end

@inline keep(c::Real, epsilon::Float64) = !iszero(c)
@inline keep(c::AbstractFloat, epsilon::Float64) = !(abs(c) <= epsilon)
@inline keep(c::Real, ctx::Context) = keep(c, ctx.epsilon)
@inline keep(c::AbstractFloat, ctx::Context) = ctx.big_epsilon === nothing ? keep(c, ctx.epsilon) : !(abs(c) <= ctx.big_epsilon)
@inline keep(c::BigFloat, ctx::Context) = !(abs(c) <= something(ctx.big_epsilon, ctx.epsilon))
@inline function finish!(a::DA, n::Int = length(a.coeffs); filter::Bool = true)
    if filter && (a.context.epsilon != 0 || a.context.big_epsilon !== nothing)
        @inbounds for i in 1:n
            keep(a.coeffs[i], a.context) || (a.coeffs[i] = zero(eltype(a.coeffs)))
        end
    end
    @inbounds while n > 1 && iszero(a.coeffs[n])
        n -= 1
    end
    a.len = n
    return a
end
@inline allocate(ctx::Context, ::Type{T}, n::Int = last(ctx.basis.ends)) where {T <: Real} = DA{T}(zeros(T, n), 1, ctx)
# Only for kernels that write every active coefficient before returning. As with
# a reused buffer, capacity beyond len is unspecified and must never be read.
@inline allocate_undef(ctx::Context, ::Type{T}, n::Int) where {T <: Real} = DA{T}(Vector{T}(undef, n), 1, ctx)
@inline function reserve!(a::DA, n::Int)
    old = length(a.coeffs)
    if n > old
        resize!(a.coeffs, n)
        @inbounds for i in (old + 1):n
            a.coeffs[i] = zero(eltype(a.coeffs))
        end
    end
    return a
end
coefftype(::DA{T}) where {T} = T
coefftype(::Type{DA{T}}) where {T} = T

DA() = DA(0.0)
DA(x::Real) = DA{typeof(float(x))}(x)
DA(a::DA) = copy(a)
function DA{T}(x::Real = zero(T)) where {T <: Real}
    ctx = ready()
    c = convert(T, x)
    return DA{T}([keep(c, ctx) ? c : zero(T)], 1, ctx)
end
function DA{T}(i::Integer, c::Real) where {T <: Real}
    ctx = ready()
    0 <= i <= ctx.basis.variables || throw(ArgumentError("Variable index out of bounds"))
    a = allocate(ctx, T, Int(i) + 1)
    a.coeffs[i + 1] = convert(T, c)
    return finish!(a)
end
DA(i::Integer, c::Real) = DA{typeof(float(c))}(i, c)
function variable(i::Integer, ::Type{T} = Float64) where {T <: Real}
    i >= 1 || throw(ArgumentError("Variable indices start at one"))
    return DA{T}(i, one(T))
end
@inline cons(a::DA) = (valid(a); @inbounds a.coeffs[1])
cons(a::Real) = a
cons(a::AbstractArray{<:DA}) = cons.(a)
cons(a::AbstractArray{<:Real}) = a
function linear(a::DA{T}) where {T}
    ctx = valid(a)
    out = zeros(T, ctx.basis.variables)
    @inbounds for i in 1:min(length(out), a.len - 1)
        out[i] = a.coeffs[i + 1]
    end
    return out
end
function linear(v::AbstractVector{<:DA})
    T = isempty(v) ? Float64 : mapreduce(coefftype, promote_type, v)
    result = zeros(T, length(v), getMaxVariables())
    for (i, p) in enumerate(v)
        result[i, :] = linear(p)
    end
    return result
end
