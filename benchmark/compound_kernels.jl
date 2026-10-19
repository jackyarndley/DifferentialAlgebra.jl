using DifferentialAlgebra, IntervalArithmetic, SpecialFunctions

include("measurement.jl")

revision = readchomp(`git -C $(pkgdir(DifferentialAlgebra)) rev-parse --short HEAD`)
println(
    (;
        julia = VERSION, revision,
        interval_arithmetic = Base.pkgversion(IntervalArithmetic),
    )
)
function report(name, f; budget = missing)
    try
        measurement = measure(f)
        println((; path = name, budget, measurement))
        return f()
    catch err
        failure = first(split(sprint(showerror, err), '\n'))
        println((; path = name, budget, failure))
    end
    return nothing
end
for budget in (0, 32 * 1024^2)
    x, = variables(1; order = 4, table_bytes = budget)
    p, q, r = 10 + 0.09x, 2 + 0.018x, 1 + 100x
    set_coefficient_tolerance!(0.01)
    for (name, f) in (
            (:integer_power, () -> p^3), (:negative_power, () -> q^(-3)),
            (:log2, () -> log2(q)), (:hypot, () -> hypot(p, q)),
            (:hypot_scalar, () -> hypot(r, 0.005)),
            (:atan_scalar, () -> atan(0.005, r)),
            (:gamma, () -> gamma(p)), (:polynomial_power, () -> p^q),
        )
        report(name, f; budget)
    end
    println(
        (;
            budget, log2_linear = coefficient(log2(q), [1]),
            cube_quadratic = coefficient(p^3, [2]),
        )
    )
end
x, y = variables(Interval{Float64}, 2; order = 4)
for (name, c) in (
        (:point, interval(3)), (:positive, interval(3, 4)),
        (:contains_zero, interval(-3, 3)),
    )
    p = name == :contains_zero ? c * x : c * x + interval(4) * y
    result = report(Symbol(name, :_norm), () -> coefficient_norm(p, 2))
    report(Symbol(name, :_degrees), () -> DifferentialAlgebra.degree_norms(p, 0, 2))
    if !isnothing(result)
        println(
            (;
                path = name, bounds = (inf(result), sup(result)),
                guaranteed = isguaranteed(result), decoration = decoration(result),
            )
        )
    end
end
for T in (Float64, BigFloat, Rational{BigInt})
    variables(T, 1; order = 2)
    A = TaylorPolynomial{T}.(T[4 0; 0 2])
    report(Symbol(:diagonal_, T), () -> DifferentialAlgebra.eigenpairs(A))
end
@assert !any(m -> nameof(m) == :GenericLinearAlgebra, values(Base.loaded_modules))
