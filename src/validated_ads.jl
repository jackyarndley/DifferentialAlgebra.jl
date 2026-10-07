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
    _degree::Int
end
domain(a::CompiledTaylorModel) = collect(deepcopy(a._coordinates.box))
remainder(a::CompiledTaylorModel) = deepcopy(a._remainder)
coefficient_type(::CompiledTaylorModel{T}) where {T} = T
nvariables(a::CompiledTaylorModel) = length(a._coordinates.box)
max_order(a::CompiledTaylorModel) = a._order
degree(a::CompiledTaylorModel) = a._degree
Base.copy(a::CompiledTaylorModel) = CompiledTaylorModel(deepcopy(a._coefficients), copy(a._exponents), deepcopy(a._remainder), deepcopy(a._coordinates), a._order, a._degree)
function Base.deepcopy_internal(a::CompiledTaylorModel, copies::IdDict)
    return get!(copies, a) do
        CompiledTaylorModel(
            Base.deepcopy_internal(a._coefficients, copies),
            Base.deepcopy_internal(a._exponents, copies),
            Base.deepcopy_internal(a._remainder, copies),
            Base.deepcopy_internal(a._coordinates, copies), a._order, a._degree
        )
    end
end
(a::CompiledTaylorModel)(point) = evaluate(a, point)

"""
    TaylorModelPatch

A Taylor-model patch with certified enclosures on a box subdomain.
`models` is a tuple of owned `CompiledTaylorModel`s;
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
nvariables(p::TaylorModelPatch) = nvariables(first(p.models))
noutputs(p::TaylorModelPatch) = length(p.models)
max_order(p::TaylorModelPatch) = max_order(first(p.models))
degree(p::TaylorModelPatch) = maximum(degree, p.models)

"""
    PiecewiseTaylorModel

A collection of Taylor-model patches on a box subdomain partition, produced by
`adaptive_map(...; estimator=IntervalBound())`. `patches` is an immutable tuple of `TaylorModelPatch`s.
`converged` means every uniform absolute error
bound met the requested tolerance. Even unresolved patches retain valid model
enclosures. Numeric point and subbox queries include remainders, reject queries
outside the physical domain and remain valid after algebra reinitialization.
Queries use a partition tree and include both closed sides of shared faces.
Refinement reevaluates the original callback and defaults to the saved estimator
and requested retained order, which can exceed `degree(map)`.
"""
struct PiecewiseTaylorModel{T <: Real, N, E <: ADSEstimator}
    _domain::Tuple{Vararg{T}}
    patches::Tuple{Vararg{TaylorModelPatch{T}}}
    nodes::Vector{N}
    _order::Int
    _scalar::Bool
    _estimator::E
    converged::Bool
end
domain(a::PiecewiseTaylorModel) = collect(deepcopy(a._domain))
nvariables(a::PiecewiseTaylorModel) = length(a._domain)
noutputs(a::PiecewiseTaylorModel) = length(first(a.patches).models)
max_order(a::PiecewiseTaylorModel) = a._order
degree(a::PiecewiseTaylorModel) = maximum(degree(m) for p in a.patches for m in p.models)
(a::PiecewiseTaylorModel)(point) = evaluate(a, point)
Base.copy(a::PiecewiseTaylorModel) = PiecewiseTaylorModel(deepcopy(a._domain), deepcopy(a.patches), deepcopy(a.nodes), a._order, a._scalar, a._estimator, a.converged)
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

This is a thin compatibility wrapper for `adaptive_map(...; estimator=IntervalBound())`.
Its preserved defaults are `order=3`, `atol=1e-6`, `splitter=:width`; the canonical
interface defaults to `order=5`, `atol=1e-8`, `splitter=:tail`.
`:width` bisects the longest side relative to the initial box, skipping fixed
coordinates. `:tail` uses nonlinear retained tails and coefficient uncertainty
from unsatisfied outputs to select an axis, falling back to relative widths.
The scalar remainder has no exact directional attribution.
`:oriented` selects a 2D exact convex polygon partition with automatic or supplied
projection directions. Direction scores are heuristic; acceptance remains rigorous.
Recompute `f` on both children: restricting a previously constructed model does
not shrink its remainder. Invalid function domains throw, including when a
loose whole-model enclosure cannot establish validity. No time integration or
verified inversion is performed. The callback must not change algebra settings.

Resource limits throw by default. `strict=false` retains valid unresolved patches
with their status and `converged=false`. Construction uses a temporary algebra
and restores the caller's algebra, including on failure. Returned snapshots own
their data and remain numerically valid independently of the global algebra.
"""
function validated_adaptive_map(
        f, box; order::Integer = 3, atol = 1.0e-6, splitter::Symbol = :width, kwargs...
    )
    return adaptive_map(f, box; estimator = IntervalBound(), order, atol, splitter, kwargs...)
end

function adaptive_map(
        f, previous::PiecewiseTaylorModel; order::Integer = max_order(previous),
        estimator::ADSEstimator = previous._estimator, kwargs...
    )
    estimator isa IntervalBound || throw(ArgumentError("Validated box refinement requires IntervalBound"))
    return interval_ads(f, domain(previous); partition = previous, order, kwargs...)
end
