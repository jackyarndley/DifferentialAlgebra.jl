# Moments of independent standardized inputs used by the uncertainty tutorials.
uniform_moment(powers) = any(isodd, powers) ? 0.0 : prod(k -> 1 / (k + 1), powers)
normal_moment(powers) = any(isodd, powers) ? 0.0 :
    prod(k -> prod(1:2:(Int(k) - 1); init = 1.0), powers)

# Work on complete coefficient lists: polynomial products would discard the
# high-degree terms needed by the covariance. Center the lists themselves so
# a lower working order or coefficient tolerance cannot alter stored moments.
function polynomial_statistics(polynomials, moment)
    terms = monomials.(polynomials)
    mean = [
        sum(m.coefficient * moment(m.exponents) for m in list; init = 0.0)
            for list in terms
    ]
    centered = [
        [
            Monomial(constant_term(p) - μ, zeros(Int, nvariables(p)));
            filter(m -> !iszero(degree(m)), list)
        ] for (p, μ, list) in zip(polynomials, mean, terms)
    ]
    covariance = [
        sum(
            a.coefficient * b.coefficient * moment(a.exponents + b.exponents)
                for a in centered[i], b in centered[j]; init = 0.0
        ) for i in eachindex(polynomials), j in eachindex(polynomials)
    ]
    return mean, covariance
end
