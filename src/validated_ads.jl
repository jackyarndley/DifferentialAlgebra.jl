# Construction and validated numeric operations live in the optional extension.
"""
    CompiledTaylorModel
    compile(model::TaylorModel)

An owned numeric snapshot of a native Taylor model. Compilation copies retained
interval coefficients, monomial exponents, remainder, domain and normalization.
`evaluate`, `enclose`, `domain` and `remainder` retain the model contract after
algebra reinitialization. It supports numeric queries, not model arithmetic or
polynomial composition. Internal arrays are read-only implementation data.
"""
struct CompiledTaylorModel{T <: Real}
    _coefficients::Vector{T}
    _exponents::Matrix{Int}
    _remainder::T
    _coordinates::ModelCoordinates{T}
    _order::Int
end
domain(a::CompiledTaylorModel) = collect(deepcopy(a._coordinates.box))
remainder(a::CompiledTaylorModel) = deepcopy(a._remainder)
coefficient_type(::CompiledTaylorModel{T}) where {T} = T
nvariables(a::CompiledTaylorModel) = length(a._coordinates.box)
max_order(a::CompiledTaylorModel) = a._order
degree(a::CompiledTaylorModel) = maximum(sum(@view(a._exponents[:, k])) for k in axes(a._exponents, 2))
Base.copy(a::CompiledTaylorModel) = CompiledTaylorModel(deepcopy(a._coefficients), copy(a._exponents), deepcopy(a._remainder), deepcopy(a._coordinates), a._order)
Base.deepcopy_internal(a::CompiledTaylorModel, copies::IdDict) = get!(() -> copy(a), copies, a)
(a::CompiledTaylorModel)(point) = evaluate(a, point)

"""
    TaylorModelPatch

A certified box-ADS leaf. `models` is a tuple of owned `CompiledTaylorModel`s;
`error_bounds` encloses each original output minus its polynomial with stored
midpoint coefficients, uniformly on this patch. It includes both the absolute
remainder and retained coefficient widths. `depth` counts bisections; `status`
is `:converged`, `:max_depth`, `:max_patches` or `:roundoff`. Data are read-only.
`domain(patch)` returns an independent copy of the physical box.
"""
struct TaylorModelPatch{T <: Real}
    models::Tuple{Vararg{CompiledTaylorModel{T}}}
    error_bounds::Tuple{Vararg{T}}
    depth::Int
    status::Symbol
end
domain(p::TaylorModelPatch) = domain(first(p.models))

"""
    PiecewiseTaylorModel

A box partition produced by `validated_adaptive_map`. `patches` is an immutable
tuple of `TaylorModelPatch`s. `converged` means every uniform absolute error
bound met the requested tolerance. Even unresolved leaves retain valid model
enclosures. Numeric point and subbox queries include remainders, reject queries
outside the physical domain and remain valid after algebra reinitialization.
"""
struct PiecewiseTaylorModel{T <: Real}
    _domain::Tuple{Vararg{T}}
    patches::Tuple{Vararg{TaylorModelPatch{T}}}
    _order::Int
    _scalar::Bool
    converged::Bool
end
domain(a::PiecewiseTaylorModel) = collect(deepcopy(a._domain))
nvariables(a::PiecewiseTaylorModel) = length(a._domain)
noutputs(a::PiecewiseTaylorModel) = length(first(a.patches).models)
max_order(a::PiecewiseTaylorModel) = a._order
degree(a::PiecewiseTaylorModel) = maximum(degree(m) for p in a.patches for m in p.models)
(a::PiecewiseTaylorModel)(point) = evaluate(a, point)
Base.copy(a::PiecewiseTaylorModel) = PiecewiseTaylorModel(deepcopy(a._domain), Tuple(TaylorModelPatch(Tuple(copy(m) for m in p.models), deepcopy(p.error_bounds), p.depth, p.status) for p in a.patches), a._order, a._scalar, a.converged)
Base.deepcopy_internal(a::PiecewiseTaylorModel, copies::IdDict) = get!(() -> copy(a), copies, a)
function Base.show(io::IO, a::PiecewiseTaylorModel)
    return print(io, "PiecewiseTaylorModel(", length(a.patches), " patches, ", noutputs(a), " outputs, ", nvariables(a), " variables, ", a.converged ? "converged" : "unresolved", ")")
end

"""
    validated_adaptive_map(
        f, box; order = 3, atol = 1.0e-6,
        max_depth = 20, max_patches = 1024, strict = true,
        names = nothing, table_bytes = 32 * 1024^2,
        splitter = :width, directions = nothing
    )

Reevaluate the original deterministic function `f(coordinates)` with native
Taylor-model inputs on each child box. `f` returns a model/real scalar or a
nonempty vector of models/reals. Domains are guaranteed decorated intervals
with Float64 or BigFloat endpoints. Each output's uniform absolute enclosure
error relative to its retained midpoint polynomial must be at most `atol`
(a positive scalar or one positive tolerance per output). This bound includes
retained coefficient rounding and the absolute remainder; samples and ordinary
DA tail estimates are not used for acceptance.

The same method is selectable through `adaptive_map(...; estimator=IntervalBound())`.
`:width` bisects the longest side relative to the initial box, skipping fixed
coordinates. `:tail` uses retained coefficient sensitivity to select an axis.
`:oriented` selects a 2D exact convex polygon partition with automatic or supplied
projection directions. Direction scores are heuristic; acceptance remains rigorous.
Recompute `f` on both children: restricting a previously constructed model does
not shrink its remainder. Invalid function domains throw, including when a
loose whole-model enclosure cannot establish validity. No time integration or
verified inversion is performed. The callback must not change algebra settings.

Resource limits throw by default. `strict=false` retains valid unresolved leaves
with their status and `converged=false`. Construction uses a temporary algebra
and restores the caller's algebra, including on failure. Returned snapshots own
their data and remain numerically valid independently of the global algebra.
"""
function validated_adaptive_map end
validated_adaptive_map(args...; kwargs...) = throw(ArgumentError("Load IntervalArithmetic to construct certified box ADS"))
