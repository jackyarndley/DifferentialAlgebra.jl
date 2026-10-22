module DifferentialAlgebraForwardDiffExt

import DifferentialAlgebra as DA
import ForwardDiff

# Only branch/domain checks inspect primal values. All surrogate arithmetic
# retains the original (possibly nested) dual numbers.
DA.continuity_primal(x::ForwardDiff.Dual) = DA.continuity_primal(ForwardDiff.value(x))
DA.taylor_solver_coefficient(::Type{<:ForwardDiff.Dual{Tag, V}}) where {Tag, V} =
    DA.taylor_solver_coefficient(V)
DA.time_exponent(p::Rational, value::ForwardDiff.Dual) =
    DA.time_exponent(p, ForwardDiff.value(value))

# Older ForwardDiff predicates inspect only the primal. A zero-valued
# coefficient with nonzero partials must remain stored; a reciprocal still
# requires an invertible primal, including through nested Duals.
DA.coefficient_iszero(c::ForwardDiff.Dual) =
    DA.coefficient_iszero(ForwardDiff.value(c)) &&
    all(DA.coefficient_iszero, ForwardDiff.partials(c))
DA.coefficient_isone(c::ForwardDiff.Dual) =
    DA.coefficient_isone(ForwardDiff.value(c)) &&
    all(DA.coefficient_iszero, ForwardDiff.partials(c))
DA.invertible_constant(c::ForwardDiff.Dual) =
    DA.invertible_constant(ForwardDiff.value(c))

# Resolve the intersections of ForwardDiff's Dual/Real and polynomial/Real
# operations without unwrapping the dual or losing its coefficient derivatives.
Base.promote_rule(
    ::Type{DA.TaylorPolynomial{T}}, ::Type{ForwardDiff.Dual{Tag, V, N}}
) where {T <: Real, Tag, V, N} =
    DA.TaylorPolynomial{promote_type(T, ForwardDiff.Dual{Tag, V, N})}
Base.promote_rule(
    ::Type{ForwardDiff.Dual{Tag, V, N}}, ::Type{DA.TaylorPolynomial{T}}
) where {T <: Real, Tag, V, N} =
    DA.TaylorPolynomial{promote_type(T, ForwardDiff.Dual{Tag, V, N})}
Base.promote_rule(
    ::Type{DA.TaylorPolynomial}, ::Type{ForwardDiff.Dual{Tag, V, N}}
) where {Tag, V, N} = DA.TaylorPolynomial{ForwardDiff.Dual{Tag, V, N}}
Base.promote_rule(
    ::Type{ForwardDiff.Dual{Tag, V, N}}, ::Type{DA.TaylorPolynomial}
) where {Tag, V, N} = DA.TaylorPolynomial{ForwardDiff.Dual{Tag, V, N}}
for op in (:+, :-, :*, :/, :(==), :<, :<=, :>, :>=, :isless, :isequal)
    @eval begin
        Base.$op(a::DA.TaylorPolynomial, b::ForwardDiff.Dual) =
            invoke(Base.$op, Tuple{DA.TaylorPolynomial, Real}, a, b)
        Base.$op(a::ForwardDiff.Dual, b::DA.TaylorPolynomial) =
            invoke(Base.$op, Tuple{Real, DA.TaylorPolynomial}, a, b)
    end
end

Base.:^(a::DA.TaylorPolynomial, b::ForwardDiff.Dual) = a^DA.scalar_argument(a, b)
Base.:^(a::ForwardDiff.Dual, b::DA.TaylorPolynomial) =
    invoke(Base.:^, Tuple{Real, DA.TaylorPolynomial}, a, b)
Base.log(a::ForwardDiff.Dual, b::DA.TaylorPolynomial) =
    invoke(Base.log, Tuple{Real, DA.TaylorPolynomial}, a, b)

# ForwardDiff defines ternary methods for each scalar family and Dual position.
# Resolve their polynomial intersections with existing unfiltered product kernels;
# adding the third operand completes the result without an intermediate cutoff.
for A in (DA.TaylorPolynomial, ForwardDiff.Dual, Real),
        B in (DA.TaylorPolynomial, ForwardDiff.Dual, Real),
        C in (DA.TaylorPolynomial, ForwardDiff.Dual, Real)
    arguments = (A, B, C)
    DA.TaylorPolynomial in arguments && ForwardDiff.Dual in arguments || continue
    multiplication = if A === DA.TaylorPolynomial
        B === DA.TaylorPolynomial ? :(DA.product(a, b; filter = false)) :
            :(DA.scaled(a, b; filter = false))
    elseif B === DA.TaylorPolynomial
        :(DA.scaled(b, a; filter = false))
    else
        :(a * b)
    end
    families = Real in arguments ? (Real, AbstractFloat, Integer, Rational, Irrational) :
        (Real,)
    for scalar in families
        left_type, right_type, addend_type = map(T -> T === Real ? scalar : T, arguments)
        @eval Base.muladd(a::$left_type, b::$right_type, c::$addend_type) =
            $multiplication + c
    end
end
# Equal tags form more specific intersections than the independent-tag methods.
Base.muladd(
    a::ForwardDiff.Dual{Tag}, b::ForwardDiff.Dual{Tag}, c::DA.TaylorPolynomial
) where {Tag} = a * b + c
Base.muladd(
    a::ForwardDiff.Dual{Tag}, b::DA.TaylorPolynomial, c::ForwardDiff.Dual{Tag}
) where {Tag} = DA.scaled(b, a; filter = false) + c
Base.muladd(
    a::DA.TaylorPolynomial, b::ForwardDiff.Dual{Tag}, c::ForwardDiff.Dual{Tag}
) where {Tag} = DA.scaled(a, b; filter = false) + c

end
