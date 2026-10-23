# Coordinate calculus and map derivatives share one coefficient transformation.
function calculus(
        a::TaylorPolynomial, v::Integer, operation::Symbol;
        order::Int = working_order(a.algebra), filter::Bool = true
    )
    ctx = valid(a); b = ctx.basis
    1 <= v <= b.variables ||
        throw(ArgumentError("Variable index $v is outside 1:$(b.variables)"))
    T = if operation === :integrate
        unit = one(coefficient_type(a))
        typeof(unit / degree_factor(unit, 1))
    else
        coefficient_type(a)
    end
    degree = b.degrees[a.len]
    output_degree = operation === :integrate ? degree + 1 : degree - 1
    capacity = b.ends[clamp(output_degree, 0, order) + 1]
    out = allocate(ctx, T, capacity)
    exponents = zeros(Int, b.variables)
    for i in 1:a.len
        c = a.coeffs[i]
        coefficient_iszero(c) && continue
        exponents .= @view b.exponents[:, i]
        degree = b.degrees[i]
        if operation === :integrate
            degree >= order && continue
            exponents[v] += 1
            out.coeffs[rank(b, exponents, degree + 1)] = c / degree_factor(c, exponents[v])
        else
            exponents[v] == 0 && continue
            power = exponents[v]; exponents[v] -= 1
            degree - 1 <= order && (out.coeffs[rank(b, exponents, degree - 1)] = degree_factor(c, power) * c)
        end
    end
    return finish!(out; filter)
end
"""
    differentiate(p, i)
    differentiate(p, counts)

Differentiate a polynomial with respect to variable `i`, or by a vector of
derivative counts. Indices start at one. For example, `differentiate(p, [2, 1])`
computes the mixed partial with two derivatives in the first variable and one
in the second. An array and an integer index differentiate each entry, preserving
the array's shape.
"""
differentiate(a::TaylorPolynomial, i::Integer) = calculus(a, i, :differentiate)
"""
    integrate(p, i)
    integrate(p, counts)

Integrate a polynomial with respect to variable `i`, choosing zero integration constant.
A vector `counts` requests repeated integrals. Terms above the working order are
discarded. The coefficient type follows the scalar division operation.
An array and an integer index integrate each entry, preserving the array's shape.
"""
integrate(a::TaylorPolynomial, i::Integer) = calculus(a, i, :integrate)
for fn in (:differentiate, :integrate)
    @eval function $fn(a::TaylorPolynomial, counts::AbstractVector{<:Integer})
        ctx = valid(a)
        remaining = checked_counts(counts, ctx.basis; limit = typemax(Int))
        out = copy(a)
        for (i, n) in enumerate(counts)
            for _ in 1:n
                remaining -= 1
                order = if $(fn === :differentiate)
                    working_order(ctx) +
                        min(remaining, ctx.basis.order - working_order(ctx))
                else
                    working_order(ctx)
                end
                out = calculus(out, i, $(QuoteNode(fn)); order, filter = false)
                iszero(out) && return out
            end
        end
        return finish!(out, min(out.len, ctx.basis.ends[working_order(ctx) + 1]))
    end
end

# Map calculus preserves the shape of ordinary arrays, including views.
differentiate(a::AbstractArray{<:TaylorPolynomial}, i::Integer) = differentiate.(a, i)
integrate(a::AbstractArray{<:TaylorPolynomial}, i::Integer) = integrate.(a, i)

"""
    gradient(p)

Return the vector of polynomial partial derivatives of `p`.
Use `constant_term(gradient(p))` for the gradient at the expansion point.
"""
gradient(a::TaylorPolynomial) = [differentiate(a, i) for i in 1:nvariables(a)]
"""
    jacobian(p)
    jacobian(map)

Return a matrix of polynomial partial derivatives.
Rows correspond to map components and columns to independent variables.
A scalar polynomial produces a one-row matrix.
"""
jacobian(a::TaylorPolynomial) = permutedims(gradient(a))
function jacobian(v::AbstractVector{<:TaylorPolynomial})
    Base.require_one_based_indexing(v)
    n = isempty(v) ? nvariables() : compatible(v).basis.variables
    return [differentiate(p, j) for p in v, j in 1:n]
end
"""
    hessian(p)
    hessian(map)

Return the matrix of second polynomial partial derivatives of `p`.
For a vector of polynomials, return a vector of Hessian matrices.
"""
function hessian(a::TaylorPolynomial)
    ctx = valid(a)
    order = min(working_order(ctx) + 1, ctx.basis.order)
    first_derivatives = [calculus(a, i, :differentiate; order, filter = false) for i in 1:nvariables(a)]
    return [differentiate(p, j) for p in first_derivatives, j in eachindex(first_derivatives)]
end
hessian(v::AbstractVector{<:TaylorPolynomial}) = hessian.(v)

"""
    DifferentialAlgebra.hessian_tensor(map)

Stack the component Hessians in an array indexed by `(variable, variable, output)`.
The map components must share one owning algebra.
"""
function hessian_tensor(v::AbstractVector{<:TaylorPolynomial})
    if isempty(v)
        T = isconcretetype(eltype(v)) ? eltype(v) : TaylorPolynomial{Float64}
        return Array{T}(undef, nvariables(), nvariables(), 0)
    end
    compatible(v)
    return stack(hessian(v); dims = 3)
end
