"""
    Monomial(coefficient, exponents)

A coefficient and its vector of nonnegative integer exponents.
`Monomial{T}` converts the coefficient to concrete type `T`. Integer exponents
are checked before conversion to `UInt32`. Retrieve monomials from a polynomial
with `monomials(p)`. Monomials own their exponent storage and compare
coefficients and exponent vectors structurally.
"""
struct Monomial{T <: Real}
    coefficient::T
    exponents::Vector{UInt32}
    function Monomial{T}(
            coefficient::T, exponents::AbstractVector{<:Integer}
        ) where {T <: Real}
        isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
        Base.require_one_based_indexing(exponents)
        all(p -> 0 <= p <= typemax(UInt32), exponents) ||
            throw(ArgumentError("Exponents must fit UInt32"))
        return new{T}(coefficient, UInt32[p for p in exponents])
    end
end
function Monomial{T}(c::Real, powers::AbstractVector{<:Integer}) where {T <: Real}
    isconcretetype(T) || throw(ArgumentError("Choose a concrete coefficient type"))
    return Monomial{T}(coefficient_convert(T, c), powers)
end
Monomial(c::T, powers::AbstractVector{<:Integer}) where {T <: Real} = Monomial{T}(c, powers)
Monomial{T}() where {T <: Real} = Monomial{T}(zero(T), zeros(UInt32, nvariables()))
Monomial() = Monomial{Float64}()
Base.:(==)(a::Monomial, b::Monomial) =
    coefficient_isequal(a.coefficient, b.coefficient) && a.exponents == b.exponents
Base.isequal(a::Monomial, b::Monomial) = a == b
Base.hash(m::Monomial, h::UInt) = hash(m.exponents, hash(m.coefficient, h))
coefficient(m::Monomial) = m.coefficient
"""
    exponents(m::Monomial)

Return a copy of the monomial exponent vector.
"""
exponents(m::Monomial) = copy(m.exponents)
"""
    degree(p)
    degree(m::Monomial)
    degree(map::CompiledMap)

Highest total degree present in a polynomial or compiled map; for a monomial,
the sum of its exponents. The zero polynomial has degree zero.
"""
degree(m::Monomial) = foldl(
    (total, power) -> Base.Checked.checked_add(total, Int(power)), m.exponents; init = 0
)
degree(p::TaylorPolynomial) = valid(p).basis.degrees[p.len]
function checked_counts(counts::AbstractVector{<:Integer}, basis::MonomialBasis; limit::Int = basis.order)
    Base.require_one_based_indexing(counts)
    length(counts) == basis.variables || throw(DimensionMismatch("Provide one exponent or count per variable"))
    total = 0
    for count in counts
        0 <= count <= limit - total || throw(ArgumentError("Exponents or counts must be nonnegative with total at most $limit"))
        total += Int(count)
    end
    return total
end
"""
    coefficient(p, exponents)
    coefficient(m::Monomial)

Return a monomial coefficient without factorial scaling. Supply one nonnegative
exponent per independent variable, with total degree at most `max_order()`.
"""
function coefficient(a::TaylorPolynomial{T}, powers::AbstractVector{<:Integer}) where {T}
    basis = valid(a).basis
    total = checked_counts(powers, basis)
    i = rank(basis, powers, total)
    return i <= a.len ? a.coeffs[i] : zero(T)
end
"""
    set_coefficient!(p, exponents, value)

Replace one monomial coefficient in `p` and return `p`. Supply exactly one
nonnegative exponent per variable. The total degree cannot exceed the
algebra's maximum order. The value is converted to the polynomial's coefficient type.
"""
function set_coefficient!(a::TaylorPolynomial{T}, jj::AbstractVector{<:Integer}, value::Real) where {T}
    ctx = valid(a)
    degree = checked_counts(jj, ctx.basis)
    i = rank(ctx.basis, jj, degree)
    reserve!(a, i)
    @inbounds for j in (a.len + 1):i
        a.coeffs[j] = zero(T)
    end
    c = coefficient_convert(T, value)
    a.coeffs[i] = keep(c, ctx) ? c : zero(T)
    return finish!(a, max(i, a.len); filter = false)
end
"""
    monomials(p)

Return nonzero monomials in increasing total degree, then basis order.
"""
function monomials(a::TaylorPolynomial)
    b = valid(a).basis
    return [
        Monomial(a.coeffs[i], @view(b.exponents[:, i]))
            for i in 1:a.len if !coefficient_iszero(a.coeffs[i])
    ]
end
"""
    monomial(p, index)
    monomial(exponents, value = 1.0)

Retrieve the `index`th nonzero monomial of `p`, or construct a polynomial with
one monomial from a vector of exponents and a coefficient.
"""
function monomial(a::TaylorPolynomial, pos::Integer)
    pos >= 1 || throw(BoundsError(a, pos))
    b = valid(a).basis
    found = 0
    for i in 1:a.len
        if !coefficient_iszero(a.coeffs[i])
            found += 1
            found == pos && return Monomial(a.coeffs[i], @view(b.exponents[:, i]))
        end
    end
    throw(BoundsError(a, pos))
end
function trim(a::TaylorPolynomial, low::Integer, high::Integer = max_order(a))
    ctx = valid(a)
    low >= 0 && high >= 0 || throw(ArgumentError("Invalid order range"))
    low > high && return zero(a)
    n = min(a.len, ctx.basis.ends[min(high, ctx.basis.order) + 1])
    out = allocate(ctx, coefficient_type(a), n)
    @inbounds for i in 1:n
        ctx.basis.degrees[i] >= low && (out.coeffs[i] = a.coeffs[i])
    end
    return finish!(out, n)
end
"""
    DifferentialAlgebra.coefficient_product(p, q)

Internal coefficient-wise product retaining the original monomials.
For example, `(1 + 2x)` and `(3 + 4x)` give `3 + 8x`. The result uses the
working order and coefficient tolerance. Both inputs must share an algebra.
"""
function coefficient_product(a::TaylorPolynomial, b::TaylorPolynomial)
    ctx = compatible(a, b)
    n = min(a.len, b.len, ctx.basis.ends[working_order(ctx) + 1])
    out = allocate(ctx, promote_type(coefficient_type(a), coefficient_type(b)), n)
    @inbounds for i in 1:n
        out.coeffs[i] = a.coeffs[i] * b.coeffs[i]
    end
    return finish!(out, n)
end
function calculus(
        a::TaylorPolynomial, v::Integer, operation::Symbol, p::Integer = 1;
        order::Int = working_order(a.algebra), filter::Bool = true
    )
    ctx = valid(a); b = ctx.basis
    1 <= v <= b.variables ||
        throw(ArgumentError("Variable index $v is outside 1:$(b.variables)"))
    0 <= p <= b.order ||
        throw(ArgumentError("Variable power $p is outside 0:$(b.order)"))
    T = if operation === :integrate
        unit = one(coefficient_type(a))
        typeof(unit / degree_factor(unit, 1))
    else
        coefficient_type(a)
    end
    degree = b.degrees[a.len]
    output_degree = operation === :integrate ? degree + 1 :
        operation === :differentiate ? degree - 1 : degree - p
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
        elseif operation === :differentiate
            exponents[v] == 0 && continue
            power = exponents[v]; exponents[v] -= 1
            degree - 1 <= order && (out.coeffs[rank(b, exponents, degree - 1)] = degree_factor(c, power) * c)
        else
            exponents[v] < p && throw(TaylorError("Polynomial is not divisible by this variable power"))
            exponents[v] -= p
            degree - p <= order && (out.coeffs[rank(b, exponents, degree - p)] = c)
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
"""
    DifferentialAlgebra.divide_variable(p, i, power = 1)

Internal variable-power division. Divide every monomial by `x[i]^power`.
Throw `TaylorError` if a retained term
is not divisible. The result uses the working order and coefficient tolerance.
"""
divide_variable(a::TaylorPolynomial, i::Integer, p::Integer = 1) = calculus(a, i, :divide, p)
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
function degree_norms(a::TaylorPolynomial, v::Integer = 0, p::Integer = 0)
    b = valid(a).basis
    0 <= v <= b.variables && p >= 0 || throw(ArgumentError("Invalid norm arguments"))
    R = typeof(abs(float(constant_term(a))))
    out = zeros(R, b.order + 1)
    @inbounds for i in 1:a.len
        degree = v == 0 ? b.degrees[i] : b.exponents[v, i]
        out[degree + 1] = max(out[degree + 1], abs(float(a.coeffs[i])))
    end
    p == 0 && return out
    # Scale each degree before taking powers, as in LinearAlgebra.norm. This
    # avoids overflowing integers and floating powers of large/small coefficients.
    sums = zeros(R, length(out))
    @inbounds for i in 1:a.len
        degree = v == 0 ? b.degrees[i] : b.exponents[v, i]
        scale = out[degree + 1]
        if isfinite(scale) && !iszero(scale)
            sums[degree + 1] += (abs(float(a.coeffs[i])) / scale)^p
        end
    end
    @inbounds for i in eachindex(out)
        isfinite(out[i]) && (out[i] *= sums[i]^(one(R) / p))
    end
    return out
end
function estimate_norms(a::TaylorPolynomial, v::Integer = 0, p::Integer = 0, n::Integer = max_order(a); errors::Bool = false)
    n >= 0 || throw(ArgumentError("Invalid estimated order"))
    norms = degree_norms(a, v, p)
    length(norms) >= 3 || throw(TaylorError("At least second order is required to estimate norms"))
    indices = findall(i -> keep(norms[i + 1], a.algebra), 1:(length(norms) - 1))
    if length(indices) < 2
        @warn "Insufficient nonzero orders to estimate coefficient norms"
        estimates = zeros(eltype(norms), n + 1)
        return errors ? (estimates, zeros(eltype(norms), min(n, length(norms) - 1) + 1)) : estimates
    end
    x = indices; y = log.(norms[indices .+ 1])
    count = length(x); sx = sum(x); sy = sum(y)
    slope = (count * sum(x .* y) - sx * sy) / (count * sum(abs2, x) - sx * sx)
    intercept = (sy - slope * sx) / count
    estimates = [exp(intercept + slope * i) for i in 0:n]
    errors || return estimates
    count = min(length(norms), n + 1)
    return estimates, max.(zero(eltype(norms)), norms[1:count] - estimates[1:count])
end
"""
    bounds(p::TaylorPolynomial)

Heuristic monomial bounds on the implicit box `[-1,1]^n`, using ordinary scalar
accumulation. Floating rounding is not enclosed, and no generating-function
truncation error is included. Wrapping its final endpoints in an interval does
not certify this calculation. Use [`enclose`](@ref) on an explicit interval box
for an outward-rounded enclosure of the stored polynomial.
"""
function bounds(a::TaylorPolynomial)
    b = valid(a).basis
    lo = hi = constant_term(a)
    for i in 2:a.len
        c = a.coeffs[i]
        if all(iseven, @view b.exponents[:, i])
            lo += min(c, zero(c)); hi += max(c, zero(c))
        else
            lo -= abs(c); hi += abs(c)
        end
    end
    return (lower = lo, upper = hi)
end
"""
    random_polynomial(T = Float64; rng = Random.default_rng(), density = 1, scaled = true)

Generate random coefficients uniformly in `[-1, 1)`. Each monomial is included
with probability `density`; when `scaled`, divide its coefficient by `2^degree`.
Pass an explicit random number generator for reproducible experiments.
"""
function random_polynomial(::Type{T} = Float64; rng = Random.default_rng(), density::Real = 1, scaled::Bool = true) where {T <: AbstractFloat}
    0 <= density <= 1 || throw(ArgumentError("Density must be between zero and one"))
    ctx = ready(); a = allocate(ctx, T)
    for i in eachindex(a.coeffs)
        if Random.rand(rng) < density
            a.coeffs[i] = 2Random.rand(rng, T) - one(T)
            scaled && (a.coeffs[i] = ldexp(a.coeffs[i], -ctx.basis.degrees[i]))
        end
    end
    return finish!(a)
end
# Map calculus preserves the shape of ordinary arrays, including views.
trim(a::AbstractArray{<:TaylorPolynomial}, low::Integer) = trim.(a, low)
trim(a::AbstractArray{<:TaylorPolynomial}, low::Integer, high::Integer) =
    trim.(a, low, high)
differentiate(a::AbstractArray{<:TaylorPolynomial}, i::Integer) = differentiate.(a, i)
integrate(a::AbstractArray{<:TaylorPolynomial}, i::Integer) = integrate.(a, i)

"""
    DifferentialAlgebra.coefficient_dot(p, q)

Internal coefficient dot product returning a scalar.
This inspection includes degrees above the working order. Both inputs must
share an algebra; coefficient tolerance is not applied to the scalar sum.
"""
function coefficient_dot(a::TaylorPolynomial, b::TaylorPolynomial)
    compatible(a, b)
    T = promote_type(coefficient_type(a), coefficient_type(b))
    result = zero(T)
    @inbounds for i in 1:min(a.len, b.len)
        result += a.coeffs[i] * b.coeffs[i]
    end
    return result
end

"""
    DifferentialAlgebra.filter_terms(p, mask)

Internal monomial-mask operation. Select stored coefficients of `p` whose
monomials occur in `mask`. The mask's
coefficient values do not scale the selected terms. The result owns its storage
and applies the coefficient tolerance while preserving stored degrees.
"""
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
"""
    DifferentialAlgebra.nterms(p)

Internal storage inspection. Count terms using the coefficient type's
exact-zero test, independently
of the working order. `size(p)` keeps Julia's scalar semantics.
"""
nterms(p::TaylorPolynomial) = (valid(p); count(!coefficient_iszero, @view p.coeffs[1:p.len]))

"""
    DifferentialAlgebra.convergence_radius(p, tolerance, power = 1)

Internal convergence heuristic. Estimate the radius at which the first omitted
degree's extrapolated coefficient
norm reaches a positive finite `tolerance`. `power` selects the degree norm used
by `estimate_norms`. This coefficient-decay heuristic is not a rigorous bound.
"""
function convergence_radius(p::TaylorPolynomial, tolerance::Real, type::Integer = 1)
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("Tolerance must be positive and finite"))
    degree = working_order(valid(p)) + 1
    estimate = estimate_norms(p, 0, type, degree)[end]
    return (tolerance / estimate)^(one(float(estimate)) / degree)
end
