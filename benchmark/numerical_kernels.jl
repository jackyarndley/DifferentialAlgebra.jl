using DifferentialAlgebra, SpecialFunctions

include("measurement.jl")

revision = readchomp(`git -C $(pkgdir(DifferentialAlgebra)) rev-parse --short HEAD`)
println((; julia = VERSION, revision, special_functions = Base.pkgversion(SpecialFunctions)))
for budget in (0, 32 * 1024^2)
    x, y, z = variables(3; order = 6, table_bytes = budget)
    p = 2 + x / 10 + y / 5 + z^2 / 20
    cubic = p^3
    for (name, f) in (
            (:exp, () -> exp(p)), (:log, () -> log(p)),
            (:quotient, () -> (1 + x) / p), (:erf, () -> erf(p)),
            (:norm, () -> coefficient_norm(p, 2)),
            (:degree_norms, () -> DifferentialAlgebra.degree_norms(p, 0, 2)),
            (:hessian, () -> hessian(cubic)),
        )
        println((; budget, path = name, measurement = measure(f)))
    end
    reciprocal_x, = variables(1; order = 3, table_bytes = budget)
    denominator = 1.0e8 + 1.0e6reciprocal_x
    set_coefficient_tolerance!(1.0e-6)
    f = () -> 1.0e8 / denominator
    measurement = measure(f)
    result = f()
    println(
        (;
            budget, path = :scaled_quotient, constant = constant_term(result),
            linear = coefficient(result, [1]), measurement,
        )
    )
    exp_x, = variables(1; order = 10, table_bytes = budget)
    exp_argument = -5 + 100exp_x
    set_coefficient_tolerance!(1.0e-4)
    f = () -> exp(exp_argument)
    measurement = measure(f)
    println((; budget, path = :filtered_exp, tenth = coefficient(f(), [10]), measurement))
    erf_x, = variables(1; order = 2, table_bytes = budget)
    erf_argument = 4 + 100erf_x
    set_coefficient_tolerance!(1.0e-4)
    f = () -> erf(erf_argument)
    measurement = measure(f)
    println((; budget, path = :filtered_erf, quadratic = coefficient(f(), [2]), measurement))
    stored_x, = variables(1; order = 4, table_bytes = budget)
    stored = stored_x^4
    with_order(2) do
        measurement = measure(() -> differentiate(stored, [2]))
        println((; budget, path = :repeated_derivative, quadratic = coefficient(differentiate(stored, [2]), [2]), measurement))
        measurement = measure(() -> hessian(stored))
        println((; budget, path = :stored_hessian, quadratic = coefficient(hessian(stored)[1, 1], [2]), measurement))
    end
end

for T in (Float32, Float64, BigFloat)
    x, y = variables(T, 2; order = 6)
    p = T(1.0e10) * (x + y)^6
    println((; type = T, path = :coefficient_norm, measurement = measure(() -> coefficient_norm(p, 2))))
    println((; type = T, path = :degree_norms, measurement = measure(() -> DifferentialAlgebra.degree_norms(p, 0, 2))))
end

x, = variables(1; order = 3)
let matrix = [2 + x one(x); one(x) 4 - x]
    println((; path = :eigenpairs, measurement = measure(() -> DifferentialAlgebra.eigenpairs(matrix))))
end
let matrix = [1 + x zero(x); zero(x) 1 - x]
    println((; path = :diagonal_eigenpairs, measurement = measure(() -> DifferentialAlgebra.eigenpairs(matrix))))
end
