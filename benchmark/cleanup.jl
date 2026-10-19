# Select the checkout with Pkg.develop in the benchmark environment.
using DifferentialAlgebra, IntervalArithmetic, LinearAlgebra
import DiffEqBase, OrdinaryDiffEqCore
include("measurement.jl")

println((; julia = VERSION, source = pkgdir(DifferentialAlgebra)))
for budget in (0, 32 * 1024^2)
    x = variables(6; order = 4, table_bytes = budget)
    p = exp(sum(x) / 6) + prod(x[1:3])
    println((; budget, path = :hessian, measurement = measure(() -> hessian(p))))
    s = TimeSeries([0.2, 1.0, zeros(19)...])
    for f in (tan, tanh)
        println((; budget, path = string(f, " TimeSeries"), measurement = measure(() -> f(s))))
        rhs(u, p, t) = [f(u[1])]
        println((; budget, path = string(f, " recorded expansion"), measurement = measure(() -> taylor_expand(rhs, [0.2], 0.0; order = 20))))
    end
    for estimator in (GuardedTail(), IntervalBound())
        build() = adaptive_map(v -> exp(v[1] + v[2]), [-0.5, -0.5], [0.5, 0.5]; estimator, order = 3, atol = 0.001, splitter = :oriented, table_bytes = budget)
        measurement = measure(build)
        fit = build()
        println((; budget, path = string(estimator, " polygon"), patches = length(fit.patches), measurement))
    end
end

const ODECore = OrdinaryDiffEqCore
const TaylorExt = Base.get_extension(DifferentialAlgebra, :DifferentialAlgebraOrdinaryDiffEqExt)
rhs(u, p, t) = [u[2], -u[1]]
series = taylor_expand(rhs, [0.0, 1.0], 0.0; order = 20)
k = [[coefficient(s, d) for s in series] for d in 0:20]
cache = TaylorExt.TaylorCache(20, rhs, nothing, zeros(2))
function measure_dense(k, cache)
    initial, out = [0.0, 1.0], zeros(2)
    dense!() = ODECore._ode_interpolant!(out, 0.5, 0.1, initial, out, k, cache, nothing, Val{0}, nothing)
    return measure(dense!)
end
println((; path = :dense_interpolation, measurement = measure_dense(k, cache)))
