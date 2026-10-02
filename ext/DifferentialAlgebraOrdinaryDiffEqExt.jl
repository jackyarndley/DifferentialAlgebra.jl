module DifferentialAlgebraOrdinaryDiffEqExt

using DifferentialAlgebra
import DiffEqBase
import OrdinaryDiffEqCore as ODECore

struct TaylorAlgorithm <: ODECore.OrdinaryDiffEqAdaptiveAlgorithm
    order::Int
end
Base.show(io::IO, alg::TaylorAlgorithm) = print(io, "TaylorMethod(", alg.order, ")")
function DifferentialAlgebra.TaylorMethod(order::Integer = 20)
    order >= 3 || throw(ArgumentError("TaylorMethod needs a time order of at least three"))
    return TaylorAlgorithm(Int(order))
end
ODECore.alg_order(alg::TaylorAlgorithm) = alg.order
ODECore.alg_adaptive_order(alg::TaylorAlgorithm) = alg.order - 2
ODECore.isfsal(::TaylorAlgorithm) = false

struct TaylorCache{F, W, U} <: ODECore.OrdinaryDiffEqConstantCache
    order::Int
    rhs::F
    workspace::W
    tmp::U
end
ODECore.get_tmp_cache(integrator, ::TaylorAlgorithm, cache::TaylorCache) = (cache.tmp,)

function ODECore.alg_cache(
        alg::TaylorAlgorithm, u, rate_prototype,
        uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2,
        f, t, dt, reltol, p, calck, ::Val{IIP}, verbose
    ) where {IIP}
    u isa Vector{<:Real} || throw(ArgumentError("TaylorMethod requires a real state vector"))
    rhs = if IIP
        (u, p, t) -> begin
            # A constant derivative may be returned alongside time series.
            du = Vector{Number}(undef, length(u))
            f(du, u, p, t)
            return du
        end
    else
        (u, p, t) -> f(u, p, t)
    end
    workspace = DifferentialAlgebra.prepare_time_expansion(rhs, u, t; order = alg.order, parameters = p)
    if coefficient_type(first(workspace.series)) <: TaylorPolynomial && !(eltype(u) <: TaylorPolynomial)
        throw(ArgumentError("Use polynomial initial states when the ODE has uncertain polynomial parameters"))
    end
    return TaylorCache(alg.order, rhs, workspace, zero.(u))
end

function store_coefficients!(k, series)
    resize!(k, max_order(first(series)) + 1)
    for degree in 0:(length(k) - 1)
        k[degree + 1] = [coefficient(s, degree) for s in series]
    end
    return nothing
end
function ODECore.initialize!(integrator, cache::TaylorCache)
    integrator.kshortsize = cache.order + 1
    store_coefficients!(integrator.k, cache.workspace.series)
    integrator.stats.nf += 2
    return nothing
end

coefficient_size(x::Real) = abs(x)
coefficient_size(x::TaylorPolynomial) = coefficient_norm(x, 1)
tolerance_at(x::Number, i) = x
tolerance_at(x::AbstractArray, i) = x[i]

function ODECore.perform_step!(integrator, cache::TaylorCache, repeat_step = false)
    (; uprev, t, dt, p, opts) = integrator
    if !repeat_step
        DifferentialAlgebra.update_time_expansion!(cache.workspace, cache.rhs, uprev, t, p)
        integrator.stats.nf += 1
    end
    u = [s(dt) for s in cache.workspace.series]
    if opts.adaptive
        error = map(eachindex(u)) do i
            scale = tolerance_at(opts.abstol, i) + tolerance_at(opts.reltol, i) *
                max(coefficient_size(uprev[i]), coefficient_size(u[i]))
            tail = max(
                coefficient_size(coefficient(cache.workspace.series[i], cache.order - 1)) * abs(dt)^(cache.order - 1),
                coefficient_size(coefficient(cache.workspace.series[i], cache.order)) * abs(dt)^cache.order
            )
            return iszero(scale) ? (iszero(tail) ? zero(tail) : oftype(tail, Inf)) : tail / scale
        end
        ODECore.set_EEst!(integrator, opts.internalnorm(error, t))
    end
    integrator.u = u
    opts.calck && store_coefficients!(integrator.k, cache.workspace.series)
    return nothing
end

# Dense output owns the coefficients of each saved step in k. It never reads
# the current cache's expansion, which may belong to a later step or callback.
function ODECore._ode_interpolant(
        Θ, dt, y0, y1, k, cache::TaylorCache,
        idxs, ::Type{Val{D}}, differential_vars
    ) where {D}
    D >= 0 || throw(ArgumentError("Derivative order must be nonnegative"))
    indices = idxs === nothing ? eachindex(y0) : idxs isa Integer ? (idxs,) : idxs
    h = Θ * dt
    result = map(indices) do i
        value = zero(y0[i])
        for degree in cache.order:-1:D
            term = k[degree + 1][i]
            for count in 0:(D - 1)
                term *= degree - count
            end
            value = muladd(value, h, term)
        end
        return value
    end
    return idxs isa Integer ? only(result) : result
end
function ODECore._ode_interpolant!(
        out, Θ, dt, y0, y1, k, cache::TaylorCache,
        idxs, deriv::Type{Val{D}}, differential_vars
    ) where {D}
    out .= ODECore._ode_interpolant(Θ, dt, y0, y1, k, cache, idxs, deriv, differential_vars)
    return out
end
function ODECore._ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::TaylorCache,
        always_calc_begin = false, allow_calc_end = true, force_calc_end = false
    )
    if length(k) < cache.order + 1 || always_calc_begin
        series = taylor_expand(cache.rhs, uprev, t; order = cache.order, parameters = p)
        store_coefficients!(k, series)
    end
    return nothing
end
DiffEqBase.interp_summary(::Type{<:TaylorCache}, dense::Bool) = dense ? "Taylor polynomial in time" : "Linear interpolation"

end
