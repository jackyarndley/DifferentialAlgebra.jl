module DifferentialAlgebraDiffEqBaseExt
using DifferentialAlgebra, DiffEqBase

DiffEqBase.value(a::DifferentialAlgebra.TaylorPolynomial) = DifferentialAlgebra.constant_term(a)

end
