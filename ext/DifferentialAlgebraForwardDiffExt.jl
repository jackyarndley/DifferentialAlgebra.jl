module DifferentialAlgebraForwardDiffExt

import DifferentialAlgebra as DA
import ForwardDiff

# Only branch/domain checks inspect primal values. All surrogate arithmetic
# retains the original (possibly nested) dual numbers.
DA.continuity_primal(x::ForwardDiff.Dual) = DA.continuity_primal(ForwardDiff.value(x))
DA.taylor_solver_coefficient(::Type{<:ForwardDiff.Dual{Tag, V}}) where {Tag, V} =
    DA.taylor_solver_coefficient(V)

# Resolve the intersections of ForwardDiff's Dual/Real and polynomial/Real
# arithmetic without unwrapping the dual or losing its coefficient derivatives.
Base.promote_rule(
    ::Type{DA.TaylorPolynomial{T}}, ::Type{ForwardDiff.Dual{Tag, V, N}}
) where {T <: Real, Tag, V, N} =
    DA.TaylorPolynomial{promote_type(T, ForwardDiff.Dual{Tag, V, N})}
Base.promote_rule(
    ::Type{ForwardDiff.Dual{Tag, V, N}}, ::Type{DA.TaylorPolynomial{T}}
) where {T <: Real, Tag, V, N} =
    DA.TaylorPolynomial{promote_type(T, ForwardDiff.Dual{Tag, V, N})}
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::DA.TaylorPolynomial, b::ForwardDiff.Dual) =
            invoke(Base.$op, Tuple{DA.TaylorPolynomial, Real}, a, b)
        Base.$op(a::ForwardDiff.Dual, b::DA.TaylorPolynomial) =
            invoke(Base.$op, Tuple{Real, DA.TaylorPolynomial}, a, b)
    end
end

end
