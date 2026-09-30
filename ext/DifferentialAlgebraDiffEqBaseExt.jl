module DifferentialAlgebraDiffEqBaseExt
using DifferentialAlgebra, DiffEqBase

DiffEqBase.value(a::DifferentialAlgebra.TaylorPolynomial) = DifferentialAlgebra.constant_term(a)
# Runge–Kutta tableaus contain scalar constants, even when the state is polynomial.
DiffEqBase.value(::Type{DifferentialAlgebra.TaylorPolynomial{T}}) where {T} = T

end
