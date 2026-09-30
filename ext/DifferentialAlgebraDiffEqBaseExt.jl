module DifferentialAlgebraDiffEqBaseExt
using DifferentialAlgebra, DiffEqBase

DiffEqBase.value(a::DifferentialAlgebra.DA) = DifferentialAlgebra.cons(a)

end
