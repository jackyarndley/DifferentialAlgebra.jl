# Shared affine substitution kernel for translation, scaling and partial evaluation.
function affine_variable(p::TaylorPolynomial, from::Integer, to::Integer, a::Real, c::Real; cutoff::Int = valid(p).basis.order)
    ctx = valid(p); b = ctx.basis
    1 <= from <= b.variables && 1 <= to <= b.variables || throw(ArgumentError("Variable index out of bounds"))
    T = promote_type(coefficient_type(p), typeof(a), typeof(c))
    result = allocate(ctx, T)
    powers = zeros(Int, b.variables)
    # Precompute scalar powers once, including 0^0 = 1 for constant terms.
    a, c = coefficient_convert(T, a), coefficient_convert(T, c)
    ap, cp = ones(T, b.order + 1), ones(T, b.order + 1)
    for k in 1:b.order
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
        choose = one(T)
        first, last = coefficient_iszero(c) ? power : 0, coefficient_iszero(a) ? 0 : power
        for k in first:last
            degree = remaining + k
            if degree <= cutoff && !coefficient_iszero(ap[k + 1]) && !coefficient_iszero(cp[power - k + 1])
                powers[to] = target + k
                result.coeffs[rank(b, powers, degree)] += coefficient * choose * ap[k + 1] * cp[power - k + 1]
            end
            k < power && (choose = choose * degree_factor(choose, power - k) / degree_factor(choose, k + 1))
        end
    end
    return finish!(result)
end
"""
    translate_variable(p, variable, a = 1, c = 0)

Replace one coordinate by `a*x[variable]+c`, retaining the other variables.
"""
translate_variable(p::TaylorPolynomial, v::Integer, a::Real = 1, c::Real = 0) = affine_variable(p, v, v, a, c)
replace_variable(p::TaylorPolynomial, from::Integer, to::Integer, a::Real = 1) = affine_variable(p, from, to, a, 0)
scale_variable(p::TaylorPolynomial, v::Integer, a::Real = 1) = affine_variable(p, v, v, a, 0)
substitute(p::TaylorPolynomial, v::Integer, value::Real = 0) = affine_variable(p, v, v, 0, value; cutoff = valid(p).cutoff)

"Dot product of corresponding polynomial coefficients, without multiplying monomials."
function coefficient_dot(a::TaylorPolynomial, b::TaylorPolynomial)
    compatible(a, b)
    T = promote_type(coefficient_type(a), coefficient_type(b))
    result = zero(T)
    @inbounds for i in 1:min(a.len, b.len)
        result += a.coeffs[i] * b.coeffs[i]
    end
    return result
end

"Retain coefficients of p whose monomials occur in mask."
function filter_terms(p::TaylorPolynomial, mask::TaylorPolynomial)
    ctx = compatible(p, mask)
    n = min(p.len, mask.len)
    result = allocate(ctx, coefficient_type(p), n)
    @inbounds for i in 1:n
        !coefficient_iszero(mask.coeffs[i]) && (result.coeffs[i] = p.coeffs[i])
    end
    return finish!(result, n)
end
function monomial(powers::AbstractVector{<:Integer}, value::Real = 1.0)
    return set_coefficient!(TaylorPolynomial(zero(value)), powers, value)
end
"Create a polynomial with the same coefficient at every monomial."
function filled(value::Real = 1.0)
    result = allocate(ready(), typeof(float(value)))
    fill!(result.coeffs, value)
    return finish!(result)
end
"Number of nonzero coefficients (Julia's size(p) retains scalar semantics)."
nterms(p::TaylorPolynomial) = (valid(p); count(!coefficient_iszero, @view p.coeffs[1:p.len]))

"Estimated radius where the first omitted order has norm tolerance; this is not a rigorous bound."
function convergence_radius(p::TaylorPolynomial, tolerance::Real, type::Integer = 1)
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("Tolerance must be positive and finite"))
    degree = valid(p).cutoff + 1
    estimate = estimate_norms(p, 0, type, degree)[end]
    return (tolerance / estimate)^(one(float(estimate)) / degree)
end
substitute(a::AbstractArray{<:TaylorPolynomial}, i::Integer, value::Real = 0) = substitute.(a, i, value)
