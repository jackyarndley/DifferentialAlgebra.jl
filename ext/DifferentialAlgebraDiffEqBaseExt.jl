module DifferentialAlgebraDiffEqBaseExt
using DifferentialAlgebra, DiffEqBase

DiffEqBase.value(a::DifferentialAlgebra.TaylorPolynomial) =
    DiffEqBase.value(DifferentialAlgebra.constant_term(a))
# Runge–Kutta tableaus contain scalar constants, even when the state is polynomial.
function scalar_value_type(::Type{T}) where {T}
    R = DiffEqBase.value(T)
    return R === T ? T : scalar_value_type(R)
end
DiffEqBase.value(::Type{DifferentialAlgebra.TaylorPolynomial{T}}) where {T} =
    scalar_value_type(T)

end
