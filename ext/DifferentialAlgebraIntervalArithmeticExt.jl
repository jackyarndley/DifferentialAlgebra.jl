module DifferentialAlgebraIntervalArithmeticExt

import DifferentialAlgebra as DA
import IntervalArithmetic as IA
import DifferentialAlgebra: TaylorPolynomial, TaylorModel
using IntervalArithmetic: Interval

# IA 1.0.12 also offers non-rigorous rounding=:none. Its configuration record
# is inspected here, without calling configure or changing process settings.
# Keep this version-specific access isolated in the optional extension.
function validated_rounding()
    IA.configuration_options.rounding === :correct ||
        throw(ArgumentError("Validated enclosure and Taylor-model arithmetic require IntervalArithmetic.configure(rounding=:correct)"))
    return nothing
end
model_context(a) = (validated_rounding(); DA.model_arithmetic(a))
compatible_models(a, b) = (validated_rounding(); DA.model_compatible(a, b))

# Skip only a valid, guaranteed exact zero. Even a thin NG/degraded zero must
# participate in arithmetic so that its metadata is not silently erased.
DA.coefficient_iszero(c::Interval) = IA.isthinzero(c) && IA.isguaranteed(c) && IA.decoration(c) == IA.com
DA.coefficient_convert(::Type{Interval{T}}, c::Union{Integer, Rational}) where {T <: IA.NumTypes} = IA.interval(T, c)
DA.coefficient_convert(::Type{Interval{T}}, c::IA.ExactReal) where {T <: IA.NumTypes} = asinterval(Interval{T}, c)
DA.coefficient_operand(::Type{Interval{T}}, c::Real) where {T <: IA.NumTypes} = DA.coefficient_convert(Interval{T}, c)
DA.degree_factor(c::Interval{T}, n::Integer) where {T <: IA.NumTypes} = IA.interval(T, n)
DA.coefficient_muladd(a::Real, b::Real, c::Interval) = DA.coefficient_convert(typeof(c), a) * DA.coefficient_convert(typeof(c), b) + c
DA.invertible_constant(c::Interval) = !IA.isnai(c) && !IA.isempty_interval(c) && !IA.in_interval(0, c)

# IA has symmetric Real promotion/arithmetic rules. Resolve their intersections
# using methods that own a DifferentialAlgebra type, never methods on IA alone.
Base.promote_rule(::Type{TaylorPolynomial{T}}, ::Type{Interval{S}}) where {T <: Real, S <: IA.NumTypes} = TaylorPolynomial{promote_type(T, Interval{S})}
Base.promote_rule(::Type{Interval{S}}, ::Type{TaylorPolynomial{T}}) where {T <: Real, S <: IA.NumTypes} = TaylorPolynomial{promote_type(T, Interval{S})}
Base.promote_rule(::Type{TaylorPolynomial}, ::Type{Interval{S}}) where {S <: IA.NumTypes} = TaylorPolynomial
Base.promote_rule(::Type{Interval{S}}, ::Type{TaylorPolynomial}) where {S <: IA.NumTypes} = TaylorPolynomial
Base.promote_rule(::Type{TaylorPolynomial{T}}, ::Type{IA.ExactReal{S}}) where {T <: Real, S <: Real} = TaylorPolynomial{promote_type(T, S)}
Base.promote_rule(::Type{TaylorPolynomial}, ::Type{IA.ExactReal{S}}) where {S <: Real} = TaylorPolynomial
TaylorPolynomial{T}(x::IA.ExactReal) where {T <: Real} = invoke(TaylorPolynomial{T}, Tuple{Real}, x)
TaylorPolynomial(x::IA.ExactReal) = TaylorPolynomial(x.value)
Base.convert(::Type{TaylorPolynomial{T}}, x::IA.ExactReal) where {T <: Real} = TaylorPolynomial{T}(x)
Base.convert(::Type{TaylorPolynomial}, x::IA.ExactReal) = TaylorPolynomial(x)
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::TaylorPolynomial{T}, b::Interval) where {T <: Real} = invoke(Base.$op, Tuple{TaylorPolynomial{T}, Real}, a, b)
        Base.$op(a::Interval, b::TaylorPolynomial) = invoke(Base.$op, Tuple{Real, TaylorPolynomial}, a, b)
    end
end

include("interval_models.jl")
include("interval_ads.jl")
include("interval_polygon_ads.jl")
include("interval_continuity.jl")

end
