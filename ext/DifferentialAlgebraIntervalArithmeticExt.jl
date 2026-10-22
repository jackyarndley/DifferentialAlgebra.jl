module DifferentialAlgebraIntervalArithmeticExt

import DifferentialAlgebra as DA
import IntervalArithmetic as IA
import DifferentialAlgebra: TaylorPolynomial, TaylorModel
using IntervalArithmetic: Interval

# IA 1.0.12 also offers non-rigorous rounding=:none. Its configuration record
# is inspected here, without calling configure or changing process settings.
# Keep this version-specific access isolated in the optional extension.
function validated_rounding()
    isdefined(IA, :configuration_options) &&
        hasproperty(IA.configuration_options, :rounding) ||
        throw(ArgumentError("Cannot verify this IntervalArithmetic rounding configuration"))
    IA.configuration_options.rounding === :correct ||
        throw(ArgumentError("Validated enclosure and Taylor-model arithmetic require IntervalArithmetic.configure(rounding=:correct)"))
    return nothing
end
model_context(a) = (validated_rounding(); DA.model_arithmetic(a))
compatible_models(a, b) = (validated_rounding(); DA.model_compatible(a, b))

function checked_interval(x::Interval; guaranteed::Bool = false)
    !IA.isnai(x) && !IA.isempty_interval(x) && IA.isbounded(x) &&
        IA.decoration(x) >= IA.def || throw(
        ArgumentError(
            "Require finite, nonempty intervals without loose-evaluation decorations"
        )
    )
    guaranteed && !IA.isguaranteed(x) && throw(
        ArgumentError(
            "Taylor models require guaranteed intervals; an NG input was supplied"
        )
    )
    return x
end
const SupportedScalar = Union{Integer, Rational, AbstractFloat, AbstractIrrational}
# Explicit interval conversion preserves decorations and guarantee flags.
asinterval(::Type{Interval{T}}, x::Interval) where {T <: IA.NumTypes} =
    IA.interval(T, x)
asinterval(::Type{Interval{T}}, x::SupportedScalar) where {T <: IA.NumTypes} =
    IA.interval(T, x)
asinterval(::Type{Interval{T}}, x::IA.ExactReal) where {T <: IA.NumTypes} =
    asinterval(Interval{T}, x.value)
asinterval(::Type{<:Interval}, x) = throw(
    ArgumentError(
        "Unsupported scalar; supply a real value or decorated Interval explicitly"
    )
)


# Skip only a valid, guaranteed exact zero. Even a thin NG/degraded zero must
# participate in arithmetic so that its metadata is not silently erased.
DA.coefficient_iszero(c::Interval) = IA.isthinzero(c) && IA.isguaranteed(c) && IA.decoration(c) == IA.com
DA.coefficient_isone(c::Interval) = IA.isthin(c) && IA.inf(c) == 1 &&
    IA.isguaranteed(c) && IA.decoration(c) == IA.com
DA.coefficient_isinteger(c::Interval) = IA.isthin(c) && isinteger(IA.inf(c)) &&
    IA.isguaranteed(c) && IA.decoration(c) == IA.com
DA.coefficient_isequal(a::Interval, b::Interval) =
    IA.isequal_interval(a, b) && IA.decoration(a) == IA.decoration(b) &&
    IA.isguaranteed(a) == IA.isguaranteed(b)
DA.coefficient_isequal(::Interval, ::Real) = false
DA.coefficient_isequal(::Real, ::Interval) = false
DA.coefficient_convert(::Type{Interval{T}}, c::Union{Integer, Rational}) where {T <: IA.NumTypes} = IA.interval(T, c)
DA.coefficient_convert(::Type{Interval{T}}, c::IA.ExactReal) where {T <: IA.NumTypes} = asinterval(Interval{T}, c)
DA.coefficient_operand(::Type{Interval{T}}, c::Real) where {T <: IA.NumTypes} = DA.coefficient_convert(Interval{T}, c)
DA.degree_factor(c::Interval{T}, n::Integer) where {T <: IA.NumTypes} = IA.interval(T, n)
DA.coefficient_muladd(a::Real, b::Real, c::Interval) = DA.coefficient_convert(typeof(c), a) * DA.coefficient_convert(typeof(c), b) + c
DA.coefficient_power(c::Interval, n::Integer) = IA.pown(c, n)
DA.taylor_solver_coefficient(::Type{<:Interval}) = false
DA.invertible_constant(c::Interval) = !IA.isnai(c) && !IA.isempty_interval(c) && !IA.in_interval(0, c)

# Normalize by a guaranteed point bound, never by an uncertain interval. The
# coefficients still participate in arithmetic, including thin degraded zeros.
@inline interval_norm_power(c, p) = p == 2 ? c * c : IA.pown(c, p)
@inline interval_norm_root(c, p) = p == 2 ? sqrt(c) : IA.rootn(c, p)

function interval_norm_values(values::AbstractVector{Interval{T}}, p::Integer) where {T}
    result = zero(Interval{T})
    if p == 0 || p == 1
        for c in values
            value = abs(c)
            result = p == 0 ? max(result, value) : result + value
        end
        return result
    end
    scale = zero(T)
    for c in values
        scale = max(scale, IA.mag(c))
    end
    scale = isfinite(scale) && !iszero(scale) ? scale : one(T)
    point = IA.interval(T, scale)
    for c in values
        # abs(c)/point is nonnegative, so its self-product is a tight square.
        result += interval_norm_power(abs(c) / point, p)
    end
    return point * interval_norm_root(result, p)
end

function interval_degree_norms(
        a::TaylorPolynomial{Interval{T}}, v::Integer, p::Integer
    ) where {T}
    basis = DA.valid(a).basis
    count = basis.order + 1
    group(i) = (v == 0 ? basis.degrees[i] : basis.exponents[v, i]) + 1
    out = zeros(Interval{T}, count)
    if p == 0 || p == 1
        for i in 1:a.len
            k = group(i)
            c = abs(a.coeffs[i])
            out[k] = p == 0 ? max(out[k], c) : out[k] + c
        end
        return out
    end
    scales = zeros(T, count)
    for i in 1:a.len
        k = group(i)
        scales[k] = max(scales[k], IA.mag(a.coeffs[i]))
    end
    # Zero and unbounded groups need no scaling. This also avoids 0/0 and Inf/Inf.
    map!(s -> isfinite(s) && !iszero(s) ? s : one(T), scales, scales)
    for i in 1:a.len
        k = group(i)
        out[k] += interval_norm_power(abs(a.coeffs[i]) / IA.interval(T, scales[k]), p)
    end
    for k in eachindex(out)
        out[k] = IA.interval(T, scales[k]) * interval_norm_root(out[k], p)
    end
    return out
end
function DA.coefficient_norm(a::TaylorPolynomial{<:Interval}, p::Real = 0)
    validated_rounding()
    DA.valid(a)
    isinteger(p) && p >= 0 ||
        throw(ArgumentError("Use a nonnegative integer coefficient norm"))
    power = p isa Integer ? p : p <= typemax(Int) ? Int(p) : BigInt(p)
    return interval_norm_values(@view(a.coeffs[1:a.len]), power)
end
function DA.degree_norms(a::TaylorPolynomial{<:Interval}, v::Integer = 0, p::Integer = 0)
    validated_rounding()
    basis = DA.valid(a).basis
    0 <= v <= basis.variables && p >= 0 ||
        throw(ArgumentError("Invalid norm arguments"))
    return interval_degree_norms(a, v, p)
end

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

include("interval_kernels.jl")
include("interval_models.jl")
include("interval_ads.jl")
include("interval_polygon_ads.jl")
include("interval_continuity.jl")

end
