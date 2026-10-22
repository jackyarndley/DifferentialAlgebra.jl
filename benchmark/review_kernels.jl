# julia --startup-file=no --project=benchmark benchmark/review_kernels.jl
using DifferentialAlgebra, SpecialFunctions, IntervalArithmetic, LinearAlgebra
using OrdinaryDiffEqCore, DiffEqBase

include("measurement.jl")

duffing(u, p, t) = [u[2], -u[1] * (1 + 0.1u[1]^2)]
energy(u) = u[2]^2 / 2 + u[1]^2 / 2 + u[1]^4 / 40
affine_map(v) = 2 + v[1] / 3 - v[2] / 5

function report(path, operation; kwargs...)
    println((; path, measurement = measure(operation; samples = 200, kwargs...)))
    return nothing
end

for budget in (0, 32 * 1024^2)
    println("\nTable budget: ", budget)
    x, = variables(1; order = 3, table_bytes = budget)
    small = 0.5 + x
    report(:small_sum, () -> small + small)
    report(:scalar_sum, () -> small + 1.0)
    report(:small_product, () -> small * small)

    for (count, order) in ((1, 20), (4, 8))
        x = variables(count; order, table_bytes = budget)
        left = 0.7 + sum(x) / count + x[1]^2 / 3
        right = 0.5 + left
        report(Symbol("hypot_", count), () -> hypot(left, right))
    end
    x, = variables(1; order = 12, table_bytes = budget)
    for center in (15, 50, 100)
        p = center + x
        report(Symbol("besselkx_", center), () -> besselkx(2, p))
    end
    x, = variables(1; order = 30, table_bytes = budget)
    p = 50 + x
    report(:besselix_50_order30, () -> besselix(2, p))

    x = taylor_models(
        fill(interval(Float64, -1 // 50, 1 // 50), 6); order = 3, table_bytes = budget
    )
    m = exp(x[1])
    report(:model_scalar_sum, () -> m + 1)
    report(:model_exp, () -> exp(x[1]))

    x, y = variables(2; order = 4, table_bytes = budget)
    initial = [1 + 0.1x, 0.1y]
    report(:polynomial_expansion, () -> taylor_expand(duffing, initial, 0.0; order = 15))
    problem = ODEProblem(duffing, initial, (0.0, 3.0))
    report(
        :polynomial_solve,
        () -> solve(problem, TaylorMethod(15); abstol = 1.0e-10, reltol = 1.0e-10)
    )
    solution = solve(problem, TaylorMethod(15); abstol = 1.0e-10, reltol = 1.0e-10)
    println(
        (;
            accepted = solution.stats.naccept,
            energy_error = coefficient_norm(energy(last(solution.u)) - energy(initial), 0),
        )
    )

    x, y = variables(2; order = 6, table_bytes = budget)
    constant = [2.0 0.3 0.1; 0.3 -1.0 0.5; 0.1 0.5 0.7]
    direction = [0.2 -0.1 0.4; -0.1 0.3 0.0; 0.4 0.0 -0.2]
    matrix = constant .+ direction .* x .+ 0.1 .* y .* (direction .+ direction')
    report(:eigenpairs, () -> DifferentialAlgebra.eigenpairs(matrix))
end

for (label, builder) in (
        (:box_build, () -> adaptive_map(affine_map, [-1, -1], [1, 1]; order = 2)),
        (
            :polygon_axes_build,
            () -> adaptive_map(
                affine_map, [-1, -1], [1, 1]; order = 2,
                splitter = :oriented, directions = :axes
            ),
        ),
        (
            :polygon_auto_build,
            () -> adaptive_map(
                affine_map, [-1, -1], [1, 1]; order = 2, splitter = :oriented
            ),
        ),
    )
    report(label, builder)
end
