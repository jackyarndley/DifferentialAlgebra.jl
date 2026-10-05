# IntervalArithmetic supplies the implementation through a package extension.
# This file only declares native storage and the dependency-free API.
struct ModelCoordinates{T <: Real}
    box::Tuple{Vararg{T}}
    center::Tuple{Vararg{T}}
    radius::Tuple{Vararg{T}}
    normalized::Tuple{Vararg{T}}
    identity::Base.RefValue{Nothing}
end

"""
    TaylorModel
    TaylorModel(value::Real, reference::TaylorModel)
    TaylorModel(p::TaylorPolynomial, absolute_remainder, reference::TaylorModel)

A native interval-coefficient polynomial, absolute interval remainder, physical
domain and affine coordinates. Load IntervalArithmetic to construct models with
[`taylor_models`](@ref). Its invariant is `f(c + r .* ξ) ∈ P(ξ) + R` throughout
the declared normalized box. The retained total degree is fixed at construction.

The three-argument constructor asserts the supplied remainder for the stored
polynomial in the reference coordinates. It cannot certify an unknown generating
function's truncation error. There is deliberately no automatic polynomial-to-model
conversion or implicit zero remainder. Constants denote their stored real values.

Construction and `copy` own coefficient storage. Accessors return independent
copies (including BigFloat endpoints); internal fields are read-only implementation
data. Reinitializing the algebra invalidates models, just as it invalidates polynomials.
Arithmetic rejects a working order different from the retained order. Evaluation
and enclosure do not depend on the working order or coefficient tolerance.
"""
struct TaylorModel{T <: Real}
    _polynomial::TaylorPolynomial{T}
    _remainder::T
    _coordinates::ModelCoordinates{T}
    _order::Int
    function TaylorModel(p::TaylorPolynomial{T}, r::T, coordinates::ModelCoordinates{T}, order::Int, ::Val{:owned}) where {T <: Real}
        return new{T}(p, r, coordinates, order)
    end
end

function model_valid(a::TaylorModel)
    ctx = valid(a._polynomial)
    degree(a._polynomial) <= a._order || throw(ArgumentError("Model polynomial exceeds its retained order"))
    return ctx
end
function model_arithmetic(a::TaylorModel)
    ctx = model_valid(a)
    ctx.cutoff == a._order || throw(ArgumentError("Taylor-model arithmetic requires working order $(a._order); use with_order"))
    return ctx
end
function model_compatible(a::TaylorModel, b::TaylorModel)
    ctx = model_arithmetic(a)
    model_arithmetic(b) === ctx || throw(ArgumentError("Models have different algebras"))
    a._order == b._order && a._coordinates.identity === b._coordinates.identity ||
        throw(ArgumentError("Models have different coordinates, normalization, domains or orders"))
    coefficient_type(a) === coefficient_type(b) || throw(ArgumentError("Model coefficient types must match"))
    return ctx
end
model_polycopy(p::TaylorPolynomial{T}) where {T} = TaylorPolynomial{T}(deepcopy(p.coeffs[1:p.len]), p.len, valid(p))

"""
    polynomial(model)

Return an independent copy of the polynomial in normalized coordinates.
This polynomial alone omits the model's absolute remainder and validity domain.
"""
polynomial(a::TaylorModel) = (model_valid(a); model_polycopy(a._polynomial))
"""
    remainder(model)

Return a copy of the absolute interval remainder. It is not a relative error or
a derivative enclosure and is not reduced by evaluating on a smaller box.
"""
remainder(a::TaylorModel) = (model_valid(a); deepcopy(a._remainder))
"""
    domain(model)

Return an independent physical validity domain: a vector of intervals for
models/box ADS or a `ConvexPolygon` for polygon ADS.
"""
domain(a::TaylorModel) = (model_valid(a); collect(deepcopy(a._coordinates.box)))
coefficient_type(::TaylorModel{T}) where {T} = T
nvariables(a::TaylorModel) = model_valid(a).basis.variables
max_order(a::TaylorModel) = (model_valid(a); a._order)
degree(a::TaylorModel) = (model_valid(a); degree(a._polynomial))
function Base.copy(a::TaylorModel)
    model_valid(a)
    c = a._coordinates
    coordinates = ModelCoordinates(deepcopy(c.box), deepcopy(c.center), deepcopy(c.radius), deepcopy(c.normalized), c.identity)
    return TaylorModel(model_polycopy(a._polynomial), deepcopy(a._remainder), coordinates, a._order, Val(:owned))
end
Base.deepcopy_internal(a::TaylorModel, copies::IdDict) = get!(() -> copy(a), copies, a)
Base.:+(a::TaylorModel) = (model_arithmetic(a); copy(a))
(a::TaylorModel)(point) = evaluate(a, point)
function Base.show(io::IO, a::TaylorModel)
    model_valid(a)
    return print(io, "TaylorModel(order = ", a._order, ", variables = ", nvariables(a), ", remainder = ", a._remainder, ")")
end

"""
    taylor_models(box; order, names = nothing, table_bytes = 32 * 1024^2)

Initialize the existing polynomial algebra and return independent physical
coordinates as Taylor models with zero remainder. `box` is a nonempty vector or
tuple of finite, nonempty, guaranteed decorated IntervalArithmetic intervals.
Only Float64 and BigFloat endpoints are supported initially. Loading the optional
dependency or extension does not initialize an algebra. As with `variables`,
calling this function invalidates previous polynomials and models.
Validated operations require IntervalArithmetic's default `rounding=:correct`.

Centers are stored binary midpoint values. Radii are outward upper bounds on
both endpoint distances. Nonfixed normalized coordinates cover `[-1,1]`; fixed
coordinates use zero radius and normalized value zero.
"""
function taylor_models end
taylor_models(args...; kwargs...) = throw(ArgumentError("Load IntervalArithmetic to construct Taylor models"))

"""
    enclose(p::TaylorPolynomial, box)
    enclose(model::TaylorModel)
    enclose(model::TaylorModel, physical_subbox)

With IntervalArithmetic loaded, enclose the **stored polynomial** on an explicit
interval box using outward-rounded monomial arithmetic and integer powers.
Ordinary floating coefficients denote their stored binary values; rational and
integer coefficients denote exact values. Existing interval flags are preserved.
No error relative to the polynomial's original generating function is inferred.

Model enclosure includes its absolute remainder. A subbox must have exactly the
physical dimension and lie inside the validity domain; normalization is outward
rounded and fixed coordinates do not divide by zero. The original remainder is
retained. `evaluate(model, physical_point)` requires one finite scalar or thin
interval per coordinate, enforces the domain and returns an interval enclosure.
"""
function enclose end
enclose(args...) = throw(ArgumentError("Load IntervalArithmetic and supply a polynomial or Taylor model with a supported domain"))

# Keep unsafe ordinary-polynomial conveniences out of the model API.
for fn in (:differentiate, :integrate, :invert)
    @eval $fn(a::TaylorModel, args...) = throw(ArgumentError($(string(fn)) * " is unsupported for TaylorModel"))
end
compile(a::TaylorModel) = throw(ArgumentError("Load IntervalArithmetic to compile TaylorModel with its remainder and domain"))
compile(a::AbstractVector{<:TaylorModel}) = throw(ArgumentError("compile is unsupported for TaylorModel"))
invert(a::AbstractVector{<:TaylorModel}) = throw(ArgumentError("Verified model inversion is unsupported"))
for op in (:(==), :<, :<=, :>, :>=, :isless, :isequal)
    @eval begin
        Base.$op(a::TaylorModel, b::TaylorModel) = throw(ArgumentError("Whole-model comparisons are unsupported"))
        Base.$op(a::TaylorModel, b::Real) = throw(ArgumentError("Whole-model comparisons are unsupported"))
        Base.$op(a::Real, b::TaylorModel) = throw(ArgumentError("Whole-model comparisons are unsupported"))
    end
end
for T in (:Float32, :Float64, :BigFloat)
    @eval Base.$T(a::TaylorModel) = throw(ArgumentError("Model-to-scalar conversion is unsupported; use evaluate"))
end
