# Select the checkout with Pkg.develop in the benchmark environment.
using DifferentialAlgebra, IntervalArithmetic, LinearAlgebra
import DiffEqBase, OrdinaryDiffEqCore
include("measurement.jl")

revision = readchomp(`git -C $(pkgdir(DifferentialAlgebra)) rev-parse --short HEAD`)
println((; julia = VERSION, revision, interval = Base.pkgversion(IntervalArithmetic), ode_core = Base.pkgversion(OrdinaryDiffEqCore)))
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

rhs(u, p, t) = [u[2], -u[1]]
problem = DiffEqBase.ODEProblem(rhs, [0.0, 1.0], (0.0, 0.1))
solution = DiffEqBase.solve(problem, TaylorMethod(20); adaptive = false, dt = 0.1)
function measure_dense(solution)
    out = zeros(2)
    return measure(() -> solution(out, 0.05))
end
println((; path = :dense_interpolation, measurement = measure_dense(solution)))
