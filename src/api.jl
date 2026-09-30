"""
    variables(n; order, table_bytes = 32 * 1024^2)
    variables(T, n; order, table_bytes = 32 * 1024^2)

Initialize an algebra of total degree `order` and return its `n` independent variables.
`T` is a concrete real coefficient type and defaults to `Float64`. The result is
an ordinary `Vector{DA{T}}`. Reinitialization invalidates existing polynomials;
use [`variable`](@ref) to retrieve a variable in the current algebra.

# Examples
```julia
x, y = variables(2; order = 6)
p = sin(x) * exp(y)
p([0.1, 0.2])
```
"""
function variables(::Type{T}, n::Integer; order::Integer, table_bytes::Integer = 32 * 1024^2) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    init(order, n; table_bytes)
    return [variable(i, T) for i in 1:n]
end
variables(n::Integer; kwargs...) = variables(Float64, n; kwargs...)

"""
    constant_term(p)

Return the constant coefficient of a polynomial, or the constant parts of an array.
For a real scalar, return the scalar itself.
"""
constant_term(p) = cons(p)

"""
    differentiate(p, i)
    differentiate(p, counts)

Differentiate a polynomial with respect to variable `i`, or by a vector of
derivative counts. Indices start at one. For example, `differentiate(p, [2, 1])`
computes the mixed partial with two derivatives in the first variable and one
in the second. An array and an integer index differentiate each entry, preserving
the array's shape. This has the same behavior as `DifferentialAlgebra.deriv`.
"""
differentiate(p, variable) = deriv(p, variable)

# Call syntax delegates to the same numeric evaluation and composition methods.
(p::DA)(point::AbstractVector{<:Real}) = evaluate(p, point)
(map::CompiledMap)(point::AbstractVector{<:Real}) = evaluate(map, point)
