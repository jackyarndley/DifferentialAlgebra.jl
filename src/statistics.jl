"All multi-indices of total degree <= order, ordered by degree then reverse lexicographically."
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
by the factorials of their exponents to recover the moments.
"""
function raw_moments(mgf::TaylorPolynomial, no::Integer)
    return lock(ALGEBRA_LOCK) do
        valid(mgf)
        0 <= no <= max_order() || throw(ArgumentError("Moment order exceeds initialized order"))
        indices = multiindices(no, nvariables())
        T = coefficient_type(mgf)
        moments = [coefficient(mgf, jj) * prod(j -> convert(T, factorial(big(j))), jj) for jj in indices]
        indices, moments
    end
end
"""
    central_moments(mgf, order)

Return multi-indices and moments after centering the moment-generating function.
The input must have constant coefficient one.
"""
function central_moments(mgf::TaylorPolynomial, no::Integer)
    return lock(ALGEBRA_LOCK) do
        valid(mgf)
        mean = linear_part(mgf)
        shift = sum(mean[i] * variable(i, coefficient_type(mgf)) for i in eachindex(mean))
        raw_moments(exp(-shift) * mgf, no)
    end
end
