# The benchmark environment selects this checkout through its [sources] entry.
using DifferentialAlgebra, IntervalArithmetic, LinearAlgebra
import DiffEqBase, OrdinaryDiffEqCore
include("measurement.jl")

for budget in (0, 32 * 1024^2)
    x = variables(6; order = 4, table_bytes = budget)
    p = exp(sum(x) / 6) + prod(x[1:3])
    println((; budget, path = :hessian, measurement = measure(() -> hessian(p))))
    s = TimeSeries([0.2, 1.0, zeros(19)...])
    for f in (tan, tanh)
        println(
            (;
                budget, path = string(f, " TimeSeries"), measurement = measure(() -> f(s)),
            )
        )
        rhs(u, p, t) = [f(u[1])]
        println(
            (;
                budget, path = string(f, " recorded expansion"),
                measurement = measure(() -> taylor_expand(rhs, [0.2], 0.0; order = 20)),
            )
        )
    end
    for estimator in (GuardedTail(), IntervalBound())
        build() = adaptive_map(
            v -> exp(v[1] + v[2]), [-0.5, -0.5], [0.5, 0.5];
            estimator, order = 3, atol = 0.001, splitter = :oriented,
            table_bytes = budget
        )
        measurement = measure(build)
        fit = build()
        println(
            (;
                budget, path = string(estimator, " polygon"),
                patches = length(fit.patches), measurement,
            )
        )
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

function inplace_growth!(du, u, p, t)
    du[1], du[2] = u[1], -u[2]
    return nothing
end
function inplace_solve()
    initial = [1.0, 2.0]
    problem = DiffEqBase.ODEProblem(inplace_growth!, initial, (0.0, 0.3))
    return DiffEqBase.solve(problem, TaylorMethod(14); adaptive = false, dt = 0.1)
end
retry_growth(u, p, t) = u
function rejected_solve()
    problem = DiffEqBase.ODEProblem(retry_growth, [1.0], (0.0, 1.0))
    return DiffEqBase.solve(
        problem, TaylorMethod(8); dt = 1.0,
        abstol = 1.0e-14, reltol = 1.0e-14
    )
end
println((; path = :inplace_solve, measurement = measure(inplace_solve)))
measurement = measure(rejected_solve)
solution = rejected_solve()
println(
    (;
        path = :rejected_solve, accepted = solution.stats.naccept,
        rejected = solution.stats.nreject, rhs_calls = solution.stats.nf, measurement,
    )
)
