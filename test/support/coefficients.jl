using DifferentialAlgebra: TaylorPolynomial, coefficient_type, nvariables

# Compare stored coefficients directly: subtraction would apply the working
# cutoff and tolerance before the test sees the difference.
coefficient_error(a::Real, b::Real) = abs(a - b)
coefficient_error(a::Union{Integer, Rational}, b::Union{Integer, Rational}) =
    abs(big(a) - big(b))
coefficientwise_error(a::Real, b::Real) = coefficient_error(a, b)
function coefficientwise_error(a::TaylorPolynomial, b::TaylorPolynomial)
    nvariables(a) == nvariables(b) || throw(DimensionMismatch("Different variable counts"))
    return maximum(1:max(a.len, b.len); init = 0) do i
        ac = i <= a.len ? a.coeffs[i] : zero(coefficient_type(a))
        bc = i <= b.len ? b.coeffs[i] : zero(coefficient_type(b))
        coefficient_error(ac, bc)
    end
end
function coefficientwise_error(a::TaylorPolynomial, b::Real)
    nvariables(a)
    return maximum(1:a.len; init = 0) do i
        coefficient_error(a.coeffs[i], i == 1 ? b : zero(b))
    end
end
coefficientwise_error(a::Real, b::TaylorPolynomial) = coefficientwise_error(b, a)
function coefficientwise_error(a::AbstractArray, b::AbstractArray)
    axes(a) == axes(b) || throw(DimensionMismatch("Different array axes"))
    return maximum(i -> coefficientwise_error(a[i], b[i]), eachindex(a, b); init = 0)
end
