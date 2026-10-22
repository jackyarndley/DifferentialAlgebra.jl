module Fixtures

using DifferentialAlgebra

export coefficients, samecoefficients, ShiftedVector, polynomial_rhs, growth!

# Inspect stored coefficients: subtracting polynomials would filter the error.
coefficients(p) = [
    coefficient(p, a) for a in
        DifferentialAlgebra.multiindices(max_order(p), nvariables(p))
]
function samecoefficients(p, q; kwargs...)
    return all(isapprox.(coefficients(p), coefficients(q); kwargs...))
end

struct ShiftedVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(v::ShiftedVector) = size(v.data)
Base.axes(v::ShiftedVector) = (0:(length(v) - 1),)
Base.getindex(::ShiftedVector, ::Int) = error("Indexed before validating axes")

function polynomial_rhs(u, p, t)
    return [
        muladd(t, t, t), muladd(t, t, 2), muladd(t, 2, t), muladd(2, t, t),
        muladd(t, 2, 2), muladd(2, t, 2), muladd(2, 2, t),
    ]
end
function growth!(du, u, p, t)
    fill!(du, 0)
    du[1] += first(p) * u[1]
    return nothing
end

end
