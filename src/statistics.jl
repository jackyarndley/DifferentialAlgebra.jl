"""
    DifferentialAlgebra.multiindices(order, nvariables)

Return vectors of `UInt32` exponents with total degree at most `order`,
including the constant term. Results are ordered by total degree, then
reverse lexicographically within each degree.
"""
function multiindices(no::Integer, nv::Integer)
    0 <= no <= 65535 && 1 <= nv <= 1024 || throw(ArgumentError("Invalid order or dimension"))
    count = binomial(big(no) + nv, nv)
    count * nv <= 32 * 1024^2 || throw(ArgumentError("Monomial basis is too large"))
    result = Vector{UInt32}[]
    sizehint!(result, Int(count))
    foreach_monomial(Int(no), Int(nv)) do powers, _
        push!(result, UInt32.(powers))
    end
    return result
end
"""
    raw_moments(mgf, order)

Return multi-indices and raw moments through total degree `order`.
`mgf` must be a moment-generating function about zero. Coefficients are multiplied
by the factorials of their exponents to recover the moments. Fixed-width integer
moments throw on overflow; use `BigInt` coefficients for larger exact moments.
"""
function raw_moments(mgf::TaylorPolynomial, no::Integer)
    basis = valid(mgf).basis
    0 <= no <= basis.order || throw(ArgumentError("Moment order exceeds initialized order"))
    indices = [UInt32.(view(basis.exponents, :, i)) for i in 1:basis.ends[no + 1]]
    moments = map(indices) do powers
        value = coefficient(mgf, powers)
        coefficient_iszero(value) && return value
        # The factorial itself may overflow even when its product with a tiny
        # coefficient is finite. Small factors also retain exact scalar arithmetic.
        for power in powers, factor in 2:Int(power)
            value = value isa Base.BitInteger ?
                Base.Checked.checked_mul(value, convert(typeof(value), factor)) :
                value * degree_factor(value, factor)
        end
        return value
    end
    return indices, moments
end
"""
    central_moments(mgf, order)

Return multi-indices and moments after centering the moment-generating function.
The input must have constant coefficient one. The requested order must not exceed
the current working order; increase it explicitly before calling this function.
This function never changes the working order.
"""
function central_moments(mgf::TaylorPolynomial, no::Integer)
    ctx = valid(mgf)
    0 <= no <= working_order(ctx) || throw(
        ArgumentError(
            "Central moment order must be between zero and the working order; " *
                "use with_order explicitly for a higher order"
        )
    )
    constant_term(mgf) == one(coefficient_type(mgf)) || throw(
        ArgumentError(
            "A moment-generating function must have constant coefficient one"
        )
    )
    mean = linear_part(mgf)
    T = negative_coefficient_type(coefficient_type(mgf))
    shift = allocate(ctx, T, length(mean) + 1)
    for i in eachindex(mean)
        shift.coeffs[i + 1] = -coefficient_convert(T, mean[i])
    end
    finish!(shift; filter = false)
    return raw_moments(product(exponential(shift; filter = false), mgf), no)
end
