# Normalized coordinates have exact powers. Retain their decorations and NG
# flags through IA operations; general boxes keep correctly rounded pown.
coordinate_power(x::Interval, k) = IA.pown(x, k)
function coordinate_power(x::Interval{Float64}, k)
    k <= 0 && return IA.pown(x, k)
    IA.isthinzero(x) && return x
    IA.inf(x) == -1 && IA.sup(x) == 1 && return isodd(k) ? x : abs(x)
    return IA.pown(x, k)
end

function finite_guaranteed_coefficients(p::TaylorPolynomial{Interval{Float64}})
    return all(@view(p.coeffs[1:p.len])) do c
        IA.isbounded(c) && !IA.isempty_interval(c) &&
            IA.decoration(c) >= IA.def && IA.isguaranteed(c)
    end
end

# Private finite, guaranteed Float64 kernel on a normalized unit box. Products
# and sums retain IA's outward rounding. Multiplication by [0,1], [-1,1] or
# {0} only selects endpoints, avoiding a second rounded interval product.
# Track every intermediate decoration: overflow can degrade a product even
# when a fixed coordinate subsequently makes that term zero.
function discarded_unit_kernel(p, q, order, bounds)
    basis = p.algebra.basis
    fixed, zero_bound, _, odd_bound = bounds
    result = zero(IA.BareInterval{Float64})
    decoration = IA.com
    @inbounds for i in 1:p.len
        ac = p.coeffs[i]
        DA.coefficient_iszero(ac) && continue
        degree = order - basis.degrees[i]
        first_j = degree < 0 ? 1 : degree >= basis.degrees[q.len] ?
            q.len + 1 : basis.ends[degree + 1] + 1
        parity, support = monomial_masks(basis, i)
        for j in first_j:q.len
            bc = q.coeffs[j]
            DA.coefficient_iszero(bc) && continue
            qparity, qsupport = monomial_masks(basis, j)
            product = IA.bareinterval(ac) * IA.bareinterval(bc)
            term = if !iszero((support | qsupport) & fixed)
                IA.bareinterval(zero_bound)
            elseif iszero(parity ⊻ qparity)
                IA.bareinterval(
                    Float64, min(0.0, IA.inf(product)),
                    max(0.0, IA.sup(product))
                )
            else
                magnitude = IA.mag(product)
                IA.bareinterval(Float64, -magnitude, magnitude)
            end
            result += term
            decoration = min(
                decoration, IA.decoration(ac), IA.decoration(bc),
                IA.decoration(odd_bound), IA.decoration(product), IA.decoration(result)
            )
        end
    end
    return IA.interval(Float64, result, decoration)
end

# Integer powers, rather than repeated interval products, preserve even powers.
# The matrix caches all powers needed by either a range or a discarded product.
function power_cache(box, order, ::Type{I}) where {I <: Interval}
    powers = Matrix{I}(undef, length(box), order + 1)
    for v in eachindex(box)
        x = asinterval(I, box[v])
        powers[v, 1] = IA.pown(x, 0)
        for k in 1:order
            powers[v, k + 1] = coordinate_power(x, k)
        end
    end
    return powers
end
function monomial_bound(coefficients, exponents, powers, ::Type{I}) where {I <: Interval}
    result = zero(I)
    # Even a constant or zero polynomial retains the domain's scalar metadata.
    for v in axes(powers, 1)
        result *= powers[v, 1]
    end
    @inbounds for k in eachindex(coefficients)
        c = asinterval(I, coefficients[k])
        checked_interval(c)
        DA.coefficient_iszero(c) && continue
        term = c
        for v in axes(exponents, 1)
            term *= powers[v, exponents[v, k] + 1]
        end
        result += term
    end
    return result
end
function polynomial_bound(p, powers, ::Type{I}) where {I <: Interval}
    b = DA.valid(p).basis
    return monomial_bound(@view(p.coeffs[1:p.len]), b.exponents, powers, I)
end

# On a normalized unit box, a nonconstant monomial ranges over [0,1] if all
# exponents are even, and [-1,1] otherwise. A fixed-axis factor makes it zero.
# Even zeroth powers carry coordinate metadata, so every axis contributes its
# decoration. Other boxes, NG coordinates and wider dimensions use powers.
function unit_product_bounds(powers, variables, ::Type{I}) where {I <: Interval}
    variables <= Sys.WORD_SIZE || return nothing
    fixed = zero(UInt)
    decoration = IA.com
    for v in 1:variables
        x = powers[v, 2]
        IA.isguaranteed(x) || return nothing
        if IA.isthinzero(x)
            fixed |= one(UInt) << (v - 1)
        elseif IA.inf(x) != -1 || IA.sup(x) != 1
            return nothing
        end
        decoration = min(decoration, IA.decoration(x))
    end
    T = IA.numtype(I)
    return (
        fixed, IA.setdecoration(zero(I), decoration),
        IA.setdecoration(IA.interval(T, 0, 1), decoration),
        IA.setdecoration(IA.interval(T, -1, 1), decoration),
    )
end

# Parent indices remove one variable at a time. Walking them uses only integer
# operations and avoids a separate per-product exponent or mask workspace.
@inline function monomial_masks(basis, k)
    parity, support = zero(UInt), zero(UInt)
    @inbounds while k > 1
        bit = one(UInt) << (basis.factors[k] - 1)
        parity ⊻= bit
        support |= bit
        k = basis.parents[k]
    end
    return parity, support
end

function discarded_unit_bound(p, q, order, bounds, ::Type{I}) where {I <: Interval}
    if I === Interval{Float64} &&
            finite_guaranteed_coefficients(p) && finite_guaranteed_coefficients(q)
        return discarded_unit_kernel(p, q, order, bounds)
    end
    return decorated_discarded_unit_bound(p, q, order, bounds, I)
end
function decorated_discarded_unit_bound(
        p, q, order, bounds, ::Type{I}
    ) where {I <: Interval}
    basis = p.algebra.basis
    fixed, zero_bound, even_bound, odd_bound = bounds
    result = zero(I)
    @inbounds for i in 1:p.len
        ac = p.coeffs[i]
        DA.coefficient_iszero(ac) && continue
        degree = order - basis.degrees[i]
        first_j = degree < 0 ? 1 : degree >= basis.degrees[q.len] ?
            q.len + 1 : basis.ends[degree + 1] + 1
        parity, support = monomial_masks(basis, i)
        for j in first_j:q.len
            bc = q.coeffs[j]
            DA.coefficient_iszero(bc) && continue
            qparity, qsupport = monomial_masks(basis, j)
            factor = if !iszero((support | qsupport) & fixed)
                zero_bound
            elseif iszero(parity ⊻ qparity)
                even_bound
            else
                odd_bound
            end
            result += (ac * bc) * factor
        end
    end
    return result
end

# Bound every ordered discarded product directly in the original basis. No
# order-2n algebra or nearest-rounded coefficient sum is used.
function discarded_bound(p, q, order, powers, ::Type{I}) where {I <: Interval}
    basis = p.algebra.basis
    basis.degrees[p.len] + basis.degrees[q.len] <= order && return zero(I)
    bounds = unit_product_bounds(powers, basis.variables, I)
    return bounds === nothing ? discarded_power_bound(p, q, order, powers, I) :
        discarded_unit_bound(p, q, order, bounds, I)
end
function discarded_power_bound(p, q, order, powers, ::Type{I}) where {I <: Interval}
    basis = p.algebra.basis
    result = zero(I)
    @inbounds for i in 1:p.len
        ac = p.coeffs[i]
        DA.coefficient_iszero(ac) && continue
        for j in 1:q.len
            basis.degrees[i] + basis.degrees[j] > order || continue
            bc = q.coeffs[j]
            DA.coefficient_iszero(bc) && continue
            term = ac * bc
            for v in 1:basis.variables
                term *= powers[v, basis.exponents[v, i] + basis.exponents[v, j] + 1]
            end
            result += term
        end
    end
    return result
end
