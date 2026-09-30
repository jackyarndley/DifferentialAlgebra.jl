struct AlgebraicVector{T<:Real} <: AbstractVector{T}
    data::Vector{T}
    AlgebraicVector{T}(data::Vector{T}) where {T<:Real} = new{T}(data)
end
function AlgebraicVector{T}(n::Integer) where {T<:Real}
    n >= 0 || throw(ArgumentError("Negative vector length"))
    AlgebraicVector{T}([zero(T) for _ in 1:n])
end
AlgebraicVector{T}(v::AbstractVector) where {T<:Real} = AlgebraicVector{T}(collect(T, v))
AlgebraicVector{T}() where {T<:Real} = AlgebraicVector{T}(T[])
function AlgebraicVector{T}(n::Integer,value::Real) where {T<:Real}
    n >= 0 || throw(ArgumentError("Negative vector length"))
    AlgebraicVector{T}([copy(convert(T,value)) for _ in 1:n])
end
AlgebraicVector(v::AbstractVector{T}) where {T<:Real} = AlgebraicVector{T}(collect(v))
Base.size(v::AlgebraicVector) = size(v.data)
Base.IndexStyle(::Type{<:AlgebraicVector}) = IndexLinear()
Base.getindex(v::AlgebraicVector, i::Int) = v.data[i]
Base.setindex!(v::AlgebraicVector, x, i::Int) = (v.data[i] = x)
Base.similar(v::AlgebraicVector, ::Type{T}, dims::Dims) where T = Array{T}(undef, dims)
Base.copy(v::AlgebraicVector) = AlgebraicVector(copy.(v.data))

struct AlgebraicMatrix{T<:Real} <: AbstractMatrix{T}
    data::Matrix{T}
    AlgebraicMatrix{T}(data::Matrix{T}) where {T<:Real} = new{T}(data)
end
AlgebraicMatrix{T}(n::Integer) where {T<:Real} = AlgebraicMatrix{T}(n, n)
AlgebraicMatrix{T}() where {T<:Real} = AlgebraicMatrix{T}(0,0)
AlgebraicMatrix{T}(m::Integer, n::Integer) where {T<:Real} = AlgebraicMatrix{T}(m, n, 0)
function AlgebraicMatrix{T}(m::Integer, n::Integer, value::Real) where {T<:Real}
    m >= 0 && n >= 0 || throw(ArgumentError("Negative matrix dimensions"))
    AlgebraicMatrix{T}([copy(convert(T, value)) for _ in 1:m, _ in 1:n])
end
AlgebraicMatrix{T}(m::AbstractMatrix) where {T<:Real} = AlgebraicMatrix{T}(Matrix{T}(m))
AlgebraicMatrix(m::AbstractMatrix{T}) where {T<:Real} = AlgebraicMatrix{T}(Matrix(m))
Base.size(m::AlgebraicMatrix) = size(m.data)
Base.IndexStyle(::Type{<:AlgebraicMatrix}) = IndexCartesian()
Base.getindex(m::AlgebraicMatrix, i::Int, j::Int) = m.data[i,j]
Base.setindex!(m::AlgebraicMatrix, x, i::Int, j::Int) = (m.data[i,j] = x)
Base.similar(m::AlgebraicMatrix, ::Type{T}, dims::Dims) where T = Array{T}(undef, dims)
Base.copy(m::AlgebraicMatrix) = AlgebraicMatrix(copy.(m.data))
Base.transpose(m::AlgebraicMatrix) = AlgebraicMatrix(permutedims(m.data))

for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::AlgebraicVector, b::AlgebraicVector) = AlgebraicVector(Base.$op.(a.data, b.data))
        Base.$op(a::AlgebraicVector, b::Real) = AlgebraicVector(Base.$op.(a.data, b))
        Base.$op(a::Real, b::AlgebraicVector) = AlgebraicVector(Base.$op.(a, b.data))
    end
end
Base.:-(a::AlgebraicVector) = AlgebraicVector(-a.data)
for f in (:sqrt,:cbrt,:exp,:log,:log2,:log10,:sin,:cos,:tan,:asin,:acos,:atan,
          :sinh,:cosh,:tanh,:asinh,:acosh,:atanh,:inv)
    @eval Base.$f(v::AlgebraicVector) = AlgebraicVector(Base.$f.(v.data))
end
sqr(v::AlgebraicVector) = AlgebraicVector(sqr.(v.data))
root(v::AlgebraicVector,p::Integer=2) = AlgebraicVector(root.(v.data,p))
isrt(v::AlgebraicVector) = root(v,-2)
icrt(v::AlgebraicVector) = root(v,-3)
Base.:^(v::AlgebraicVector,p::Real) = AlgebraicVector(v.data.^p)
Base.atan(y::AlgebraicVector,x::AlgebraicVector) = AlgebraicVector(atan.(y.data,x.data))
Base.log(b::Real,v::AlgebraicVector) = AlgebraicVector(log.(b,v.data))
Base.log(::Irrational{:ℯ},v::AlgebraicVector) = log(v)
vnorm(v::AbstractVector{<:Real}) = sqrt(sum(abs2,v))
LinearAlgebra.normalize(v::AlgebraicVector,p::Real=2) = p == 2 ? v/vnorm(v) : throw(ArgumentError("Use the Euclidean norm"))
extract(v::AlgebraicVector,first::Integer,last::Integer) = AlgebraicVector(copy.(v.data[first:last]))
concat(v::AlgebraicVector,w::AbstractVector{<:Real}) = AlgebraicVector(vcat(v.data,w))
cons(a::AbstractArray{<:DA}) = cons.(a)
cons(a::AbstractArray{<:Real}) = a
trim(v::AbstractVector{<:DA}, low::Integer, high::Integer=getMaxOrder()) = trim.(v, low, high)
deriv(v::AbstractVector{<:DA}, i::Integer) = deriv.(v, i)
integrate(v::AbstractVector{<:DA}, i::Integer) = integrate.(v, i)
plug(v::AbstractVector{<:DA},i::Integer,value::Real=0) = plug.(v,i,value)
toString(a::Union{AlgebraicVector,AlgebraicMatrix}) = sprint(show, MIME"text/plain"(), a.data)

"Return the first n independent variables."
function identity(n::Integer=getMaxVariables())
    0 <= n <= getMaxVariables() || throw(ArgumentError("Invalid identity dimension"))
    AlgebraicVector(variable.(1:n))
end
function identity(::Type{T},n::Integer=getMaxVariables()) where {T<:Real}
    0 <= n <= getMaxVariables() || throw(ArgumentError("Invalid identity dimension"))
    AlgebraicVector([variable(i,T) for i in 1:n])
end
function identity(indices::AbstractVector{<:Integer}, sorted::Bool=false)
    all(i -> 1 <= i <= getMaxVariables(), indices) || throw(ArgumentError("Variable index out of bounds"))
    sorted || return AlgebraicVector(variable.(indices))
    result = AlgebraicVector{DA}(getMaxVariables())
    for i in indices
        result[i] = variable(i)
    end
    result
end
