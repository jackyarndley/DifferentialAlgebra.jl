module DifferentialAlgebraOrdinaryDiffEqExt

using DifferentialAlgebra
import DiffEqBase
import OrdinaryDiffEqCore as ODECore

struct TaylorAlgorithm <: ODECore.OrdinaryDiffEqAdaptiveAlgorithm
    order::Int
end
Base.show(io::IO, alg::TaylorAlgorithm) = print(io, "TaylorMethod(", alg.order, ")")
function DifferentialAlgebra.TaylorMethod(order::Integer)
    order >= 3 || throw(ArgumentError("TaylorMethod needs a time order of at least three"))
    return TaylorAlgorithm(Int(order))
end
ODECore.alg_order(alg::TaylorAlgorithm) = alg.order
ODECore.alg_adaptive_order(alg::TaylorAlgorithm) = alg.order - 2
ODECore.isfsal(::TaylorAlgorithm) = false

# SciML dispatches on distinct cache supertypes for allocating and in-place
# problems. Keep those concrete containers; their stepping/interpolation logic
# is shared through TaylorCaches below.
mutable struct TaylorCache{F, W, U, T} <: ODECore.OrdinaryDiffEqConstantCache
    order::Int
    rhs::F
    workspace::W
    tmp::U
    expansion_step::Int
    expansion_time::T
end
mutable struct TaylorMutableCache{F, W, U, T} <: ODECore.OrdinaryDiffEqMutableCache
    order::Int
    rhs::F
    workspace::W
    tmp::U
    expansion_step::Int
    expansion_time::T
end
const TaylorCaches = Union{TaylorCache, TaylorMutableCache}
ODECore.full_cache(cache::TaylorMutableCache) = (cache.tmp,)
ODECore.get_fsalfirstlast(::TaylorMutableCache, u) = (zero(u), zero(u))
ODECore.get_tmp_cache(integrator, ::TaylorAlgorithm, cache::TaylorCache) = (cache.tmp,)

# Validate common parameter containers before an in-place RHS stores a polynomial
# or unsupported interval coefficient in a numeric derivative buffer.
check_parameters(p, polynomial_state) = nothing
function check_parameters(p::Real, polynomial_state)
    DifferentialAlgebra.taylor_solver_coefficient(typeof(p)) ||
        throw(ArgumentError("TaylorMethod does not support interval parameters"))
    p isa TaylorPolynomial && !polynomial_state && throw(
        ArgumentError(
            "Use polynomial initial states when the ODE has uncertain polynomial parameters"
        )
    )
    return nothing
end
function check_parameters(p::Union{AbstractArray, Tuple, NamedTuple}, polynomial_state)
    foreach(value -> check_parameters(value, polynomial_state), p)
    return nothing
end

function ODECore.alg_cache(
        alg::TaylorAlgorithm, u, rate_prototype,
        uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2,
        f, t, dt, reltol, p, calck, ::Val{IIP}, verbose
    ) where {IIP}
    u isa Vector{<:Real} || throw(ArgumentError("TaylorMethod requires a real state vector"))
    all(value -> DifferentialAlgebra.taylor_solver_coefficient(typeof(value)), u) ||
        throw(ArgumentError("TaylorMethod does not support interval coefficients"))
    check_parameters(p, eltype(u) <: TaylorPolynomial)
    rhs = if IIP
        (u, p, t) -> begin
            du = DifferentialAlgebra.time_rhs_buffer(u)
            f(du, u, p, t)
            return du
        end
    else
        (u, p, t) -> f(u, p, t)
    end
    workspace = DifferentialAlgebra.prepare_time_expansion(rhs, u, t; order = alg.order, parameters = p)
    DifferentialAlgebra.taylor_solver_coefficient(coefficient_type(first(workspace.series))) ||
        throw(ArgumentError("TaylorMethod does not support interval coefficients"))
    if coefficient_type(first(workspace.series)) <: TaylorPolynomial && !(eltype(u) <: TaylorPolynomial)
        throw(ArgumentError("Use polynomial initial states when the ODE has uncertain polynomial parameters"))
    end
    return IIP ? TaylorMutableCache(alg.order, rhs, workspace, zero.(u), -1, t) :
        TaylorCache(alg.order, rhs, workspace, zero.(u), -1, t)
end

function store_coefficients!(k, series)
    resize!(k, max_order(first(series)) + 1)
    for degree in 0:(length(k) - 1)
        k[degree + 1] = [coefficient(s, degree) for s in series]
    end
    return nothing
end
function ODECore.initialize!(integrator, cache::TaylorCaches)
    cache.expansion_step = -1
    integrator.kshortsize = cache.order + 1
    store_coefficients!(integrator.k, cache.workspace.series)
    integrator.stats.nf += 2
    return nothing
end

coefficient_size(x::Real) = abs(x)
coefficient_size(x::TaylorPolynomial) = coefficient_norm(x, 1)
tolerance_at(x::Number, i) = x
tolerance_at(x::AbstractArray, i) = x[i]
function tail_size(coefficient, h, degree)
    value = coefficient_size(coefficient)
    # Starting at the coefficient avoids overflow in h^degree when their
    # product is representable, as in a tiny coefficient over a long interval.
    for _ in 1:degree
        value *= h
    end
    return value
end

function ODECore.perform_step!(integrator, cache::TaylorCaches, _repeat_step = false)
    (; uprev, t, dt, p, opts) = integrator
    # The expansion is independent of dt. OrdinaryDiffEq retries without passing
    # repeat_step, so key the cache to the accepted-step counter and start time.
    # Public state/parameter modifications invalidate it through the solver flags.
    if integrator.iter == 1 || cache.expansion_step != integrator.success_iter ||
            cache.expansion_time != t || integrator.reeval_fsal ||
            integrator.derivative_discontinuity
        DifferentialAlgebra.update_time_expansion!(cache.workspace, cache.rhs, uprev, t, p)
        integrator.stats.nf += 1
        cache.expansion_step = integrator.success_iter
        cache.expansion_time = t
    end
    u = [s(dt) for s in cache.workspace.series]
    if opts.adaptive
        error = map(eachindex(u)) do i
            scale = tolerance_at(opts.abstol, i) + tolerance_at(opts.reltol, i) *
                max(coefficient_size(uprev[i]), coefficient_size(u[i]))
            tail = max(
                tail_size(coefficient(cache.workspace.series[i], cache.order - 1), abs(dt), cache.order - 1),
                tail_size(coefficient(cache.workspace.series[i], cache.order), abs(dt), cache.order)
            )
            return iszero(scale) ? (iszero(tail) ? zero(tail) : oftype(tail, Inf)) : tail / scale
        end
        ODECore.set_EEst!(integrator, opts.internalnorm(error, t))
    end
    if cache isa TaylorMutableCache
        copyto!(integrator.u, u)
    else
        integrator.u = u
    end
    opts.calck && store_coefficients!(integrator.k, cache.workspace.series)
    return nothing
end

# Dense output owns the coefficients of each saved step in k. It never reads
# the current cache's expansion, which may belong to a later step or callback.
function dense_value(k, order, i, h, prototype, ::Val{D}) where {D}
    value = zero(prototype)
    for degree in order:-1:D
        term = k[degree + 1][i]
        for count in 0:(D - 1)
            term *= degree - count
        end
        value = muladd(value, h, term)
    end
    return value
end
function ODECore._ode_interpolant(
        Θ, dt, y0, y1, k, cache::TaylorCaches,
        idxs, ::Type{Val{D}}, differential_vars
    ) where {D}
    D >= 0 || throw(ArgumentError("Derivative order must be nonnegative"))
    h = Θ * dt
    y0 isa Number && idxs === nothing &&
        return dense_value(k, cache.order, 1, h, y0, Val(D))
    indices = idxs === nothing ? eachindex(y0) : idxs isa Integer ? (idxs,) : idxs
    idxs isa Integer && return dense_value(k, cache.order, idxs, h, y0[idxs], Val(D))
    return map(i -> dense_value(k, cache.order, i, h, y0[i], Val(D)), indices)
end
function ODECore._ode_interpolant!(
        out, Θ, dt, y0, y1, k, cache::TaylorCaches,
        idxs, deriv::Type{Val{D}}, differential_vars
    ) where {D}
    D >= 0 || throw(ArgumentError("Derivative order must be nonnegative"))
    indices = idxs === nothing ? eachindex(y0) : idxs isa Integer ? (idxs,) : idxs
    length(out) == length(indices) || throw(DimensionMismatch("Output length must match the interpolation selection"))
    h = Θ * dt
    for (j, i) in zip(eachindex(out), indices)
        out[j] = dense_value(k, cache.order, i, h, y0[i], Val(D))
    end
    return out
end
function ODECore._ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::TaylorCaches,
        always_calc_begin = false, allow_calc_end = true, force_calc_end = false
    )
    if length(k) < cache.order + 1 || always_calc_begin
        series = taylor_expand(cache.rhs, uprev, t; order = cache.order, parameters = p)
        store_coefficients!(k, series)
    end
    return nothing
end
DiffEqBase.interp_summary(::Type{<:TaylorCaches}, dense::Bool) =
    dense ? "Taylor polynomial in time" : "Linear interpolation"

end
