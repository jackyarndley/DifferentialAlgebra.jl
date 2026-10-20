module DifferentialAlgebraForwardDiffExt

import DifferentialAlgebra as DA
import ForwardDiff

# Only branch/domain checks inspect primal values. All surrogate arithmetic
# retains the original (possibly nested) dual numbers.
DA.continuity_primal(x::ForwardDiff.Dual) = DA.continuity_primal(ForwardDiff.value(x))
DA.taylor_solver_coefficient(::Type{<:ForwardDiff.Dual{Tag, V}}) where {Tag, V} =
    DA.taylor_solver_coefficient(V)

end
