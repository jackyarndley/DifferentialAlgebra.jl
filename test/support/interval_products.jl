# Exact rational interval arithmetic, independent of IA and the polynomial
# basis, supplies ordered-product remainder bounds without floating rounding.
function rational_interval_product(a, b)
    values = (a[1] * b[1], a[1] * b[2], a[2] * b[1], a[2] * b[2])
    return extrema(values)
end

function rational_discarded_product(terms_a, terms_b, box, order)
    lower, upper = big(0) // 1, big(0) // 1
    for (alpha, ac) in terms_a, (beta, bc) in terms_b
        exponents = alpha .+ beta
        sum(exponents) > order || continue
        term = rational_interval_product(big.(ac), big.(bc))
        for (exponent, endpoints) in zip(exponents, box)
            lo, hi = big.(endpoints)
            power = if iszero(exponent)
                (one(lo), one(hi))
            elseif iseven(exponent) && lo <= 0 <= hi
                (zero(lo), max(lo^exponent, hi^exponent))
            else
                extrema((lo^exponent, hi^exponent))
            end
            term = rational_interval_product(term, power)
        end
        lower += term[1]
        upper += term[2]
    end
    return lower, upper
end
