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
