# Shared affine substitution kernel for translation, scaling and partial evaluation.
function affine_variable(p::TaylorPolynomial, from::Integer, to::Integer, a::Real, c::Real; cutoff::Int = valid(p).basis.order)
    ctx = valid(p); b = ctx.basis
    1 <= from <= b.variables && 1 <= to <= b.variables || throw(ArgumentError("Variable index out of bounds"))
    T = promote_type(coefficient_type(p), typeof(a), typeof(c))
    degree = b.degrees[p.len]
    result = allocate(ctx, T, b.ends[min(degree, cutoff) + 1])
    powers = zeros(Int, b.variables)
    # Precompute scalar powers once, including 0^0 = 1 for constant terms.
    a, c = coefficient_convert(T, a), coefficient_convert(T, c)
    ap, cp = ones(T, degree + 1), ones(T, degree + 1)
    for k in 1:degree
        ap[k + 1], cp[k + 1] = ap[k] * a, cp[k] * c
    end
    for i in 1:p.len
        coefficient = coefficient_operand(T, p.coeffs[i])
        coefficient_iszero(coefficient) && continue
        powers .= @view b.exponents[:, i]
        power = powers[from]
        powers[from] = 0
        target = powers[to]
        remaining = b.degrees[i] - power
        choose = T <: Union{Integer, Rational} ? big(1) : one(T)
        first, last = coefficient_iszero(c) ? power : 0, coefficient_iszero(a) ? 0 : power
        for k in first:last
            degree = remaining + k
            if degree <= cutoff && !coefficient_iszero(ap[k + 1]) && !coefficient_iszero(cp[power - k + 1])
                powers[to] = target + k
                result.coeffs[rank(b, powers, degree)] += coefficient * choose * ap[k + 1] * cp[power - k + 1]
            end
            if k < power
                choose = if T <: Union{Integer, Rational}
                    div(choose * (power - k), k + 1)
                else
                    choose * degree_factor(choose, power - k) / degree_factor(choose, k + 1)
                end
            end
        end
    end
    return finish!(result)
end
"""
    translate_variable(p, variable, a = 1, c = 0)

Replace one coordinate by `a*x[variable]+c`, retaining the other variables.
Binomial factors are exact for integer and rational coefficients. Scalar powers
still use the promoted coefficient type: fixed-width integer powers can overflow
as in ordinary Julia arithmetic. Use `BigInt` or `Rational{BigInt}` coefficients
when unbounded exact affine arithmetic is required.
"""
translate_variable(p::TaylorPolynomial, v::Integer, a::Real = 1, c::Real = 0) = affine_variable(p, v, v, a, c)
"""
    DifferentialAlgebra.replace_variable(p, from, to, scale = 1)

Internal affine-substitution helper. Replace `x[from]` by `scale*x[to]`,
retaining the other coordinates. This uses
the affine-substitution coefficient rules of `translate_variable`.
"""
replace_variable(p::TaylorPolynomial, from::Integer, to::Integer, a::Real = 1) = affine_variable(p, from, to, a, 0)
"""
    DifferentialAlgebra.scale_variable(p, i, scale = 1)

Internal affine-substitution helper. Replace `x[i]` by `scale*x[i]`, retaining
the other coordinates. This uses the
affine-substitution coefficient rules of `translate_variable`.
"""
scale_variable(p::TaylorPolynomial, v::Integer, a::Real = 1) = affine_variable(p, v, v, a, 0)
substitute(p::TaylorPolynomial, v::Integer, value::Real = 0) =
    affine_variable(p, v, v, 0, value; cutoff = working_order(valid(p)))

substitute(a::AbstractArray{<:TaylorPolynomial}, i::Integer, value::Real = 0) = substitute.(a, i, value)
