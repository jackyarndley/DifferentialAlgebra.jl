using Test, DifferentialAlgebra, OrdinaryDiffEqCore, SciMLBase, IntervalArithmetic
using ForwardDiff, DiffEqBase
using ..Fixtures

@testset "Time-series ownership and exact recorded operations" begin
    values = [1.0, 2.0, 3.0]
    series = TimeSeries(view(values, 1:2))
    values[1] = 9
    @test series.coefficients == [1, 2]
    @test TimeSeries(1:3).coefficients == [1, 2, 3]
    expected = (
        [0, 0, 1 // 2, 1 // 3, 0], [0, 2, 0, 1 // 3, 0],
        [0, 0, 3 // 2, 0, 0], [0, 0, 3 // 2, 0, 0],
        [0, 2, 1, 0, 0], [0, 2, 1, 0, 0], [0, 4, 1 // 2, 0, 0],
    )
    for T in (Float64, Rational{BigInt})
        solution = taylor_expand(polynomial_rhs, zeros(T, 7), zero(T); order = 4)
        @test all(s.coefficients == T.(c) for (s, c) in zip(solution, expected))
    end
    solution = only(taylor_expand((u, p, t) -> u .^ 2, [1.0], 0.0; order = 8))
    @test solution.coefficients == ones(9) # u'=u², u(0)=1 -> 1/(1-t).
    x, = variables(1; order = 2)
    solution = only(
        taylor_expand(
            (u, p, t) -> p .* u, [one(x)], 0.0;
            order = 6, parameters = 1 + x
        )
    )
    @test hcat(coefficients.(solution.coefficients)...) ≈
        [binomial(k, i) / factorial(k) for i in 0:2, k in 0:6]
    @test_throws MethodError TaylorMethod(20.0)
    @test DiffEqBase.value(TaylorPolynomial{typeof(ForwardDiff.Dual(2.0, 1.0))}) === Float64
    # Closed-form coefficients check eager and recorded paths independently.
    # Non-unit expansion points exercise reciprocal factors in log and powers.
    clock = TimeSeries([0.0, 1.0, 0.0, 0.0, 0.0])
    analytic_cases = (
        (t -> exp(2t), [1, 2, 2, 4 / 3, 2 / 3]),
        (t -> log(2 + t), [log(2), 1 / 2, -1 / 8, 1 / 24, -1 / 64]),
        (t -> sqrt(4 + t), [2, 1 / 4, -1 / 64, 1 / 512, -5 / 16384]),
        (tan, [0, 1, 0, 1 / 3, 0]),
        (tanh, [0, 1, 0, -1 / 3, 0]),
        (t -> (1 + t)^(3 // 2), [1, 3 / 2, 3 / 8, -1 / 16, 3 / 128]),
        (t -> 1 / (2 - t), [1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32]),
    )
    for (f, expected) in analytic_cases
        @test f(clock).coefficients ≈ expected
        integral = only(taylor_expand((u, p, t) -> [f(t)], [0.0], 0.0; order = 5))
        @test integral.coefficients ≈ [0; expected ./ (1:5)]
    end
    shifted = only(taylor_expand((u, p, t) -> [t^2], [8.0], 2.0; order = 3))
    @test shifted.coefficients ≈ [8, 4, 2, 1 / 3]
    @test shifted(0.5) ≈ 8 + (2.5^3 - 2^3) / 3 # Evaluate a local offset.

    @test_throws ArgumentError TimeSeries(Float64[])
    @test_throws ArgumentError TimeSeries([1.0]) + clock
    @test_throws DimensionMismatch taylor_expand((u, p, t) -> [1, 2], [0.0], 0.0)
    @test_throws DifferentialAlgebra.TaylorError taylor_expand((u, p, t) -> [sqrt(t)], [0.0], 0.0)
    @test iszero(sqrt(TimeSeries([0.0])))
    # Nonzero partials do not make a zero primal invertible, even for nested duals.
    tangent = ForwardDiff.Dual(0.0, 1.0)
    for center in (tangent, ForwardDiff.Dual{:outer}(tangent, zero(tangent)))
        @test iszero(sqrt(TimeSeries([zero(center)])))
        series = TimeSeries([center, one(center)])
        for f in (inv, log, sqrt)
            @test_throws DifferentialAlgebra.TaylorError f(series)
        end
        @test_throws DifferentialAlgebra.TaylorError series^(3 // 2)
        @test_throws DifferentialAlgebra.TaylorError sqrt(TimeSeries([center]))
        @test_throws DifferentialAlgebra.TaylorError taylor_expand(
            (u, p, t) -> [sqrt(t)], [zero(center)], center; order = 1
        )
    end
end

@testset "Time recurrences retain finite scales and coefficient perturbations" begin
    # Check every coefficient separately: a large center must not hide a lost
    # small derivative in a vector norm. Recorded results integrate these RHSs.
    cases = (
        (t -> log(1.0e308 * (1 + t)), [log(1.0e308), 1, -1 / 2, 1 / 3]),
        (t -> (1.0e308 * (1 + t))^0.25, 1.0e77 .* [1, 1 / 4, -3 / 32]),
        (t -> (1.0e308 + 1.0e-100t)^0.75, [1.0e231, 7.5e-178]),
        (t -> exp(1.0e308t^2), [1, 0, 1.0e308]),
        (t -> sin(1.0e308t^2), [0, 0, 1.0e308]),
        (t -> 1.0e-320 / (1.0e-320 * (1 + t)), [1, -1, 1]),
        (t -> 1.0e300 / (1.0e200 + 1.0e300t), [1.0e100, -1.0e200, 1.0e300]),
        (t -> (-2.0)^zero(t), [1, 0, 0]),
    )
    for (f, expected) in cases
        clock = TimeSeries([0.0, 1.0, zeros(length(expected) - 2)...])
        @test all(isapprox.(f(clock).coefficients, expected; rtol = 2.0e-14, atol = 0))
        integral = only(taylor_expand((u, p, t) -> [f(t)], [0.0], 0.0; order = length(expected)))
        @test all(isapprox.(integral.coefficients, [0; expected ./ (1:length(expected))]; rtol = 2.0e-14, atol = 0))
    end
    for f in (sqrt, cbrt, t -> t^0.25)
        @test f(TimeSeries(zeros(3))).coefficients == zeros(3)
        @test iszero(only(taylor_expand((u, p, t) -> [f(zero(t))], [0.0], 0.0; order = 3)))
    end
    @test ForwardDiff.derivative(p -> coefficient(TimeSeries([2.0, 1.0])^p, 1), 2.0) ≈ 2 + 4log(2)
    @test ForwardDiff.derivative(
        p -> coefficient(only(taylor_expand((u, p, t) -> [(2 + t)^p], [0.0], 0.0; order = 2, parameters = p)), 2), 2.0
    ) ≈ 1 + 2log(2)
    extreme = coefficient(TimeSeries([1.0, 0.0, 1.0e-308])^ForwardDiff.Dual(1.0e308, 1.0), 2)
    @test ForwardDiff.value(extreme) ≈ 1 && ForwardDiff.partials(extreme)[1] ≈ 1.0e-308
    setprecision(256) do
        expected = sqrt(BigFloat(1 // 3))
        @test coefficient((1 // 3)^TimeSeries([big"0.5", big"0"]), 0) == expected
        @test coefficient(only(taylor_expand((u, p, t) -> [(1 // 3)^(zero(t) + big"0.5")], [big"0"], big"0"; order = 2)), 1) == expected
    end
    for denominator in (0.0, ForwardDiff.Dual(0.0, 1.0))
        @test_throws DifferentialAlgebra.TaylorError TimeSeries([1.0, 2.0]) / denominator
    end
    @test_throws ArgumentError taylor_expand(
        (u, p, t) -> t isa Real ? [zero(t)] : ShiftedVector([zero(t)]), [0.0], 0.0; order = 1
    )

    x, = variables(1; order = 2)
    for wrap in (identity, v -> v + zero(x))
        a = TimeSeries(wrap.([1.0e-320, 1.0, 0.0]))
        @test constant_term.((a / a).coefficients) == [1, 0, 0]
        a = TimeSeries(wrap.([1.0e308, 0.0, 1.0e308]))
        b = TimeSeries(wrap.([1.0e308, -1.0e308, 1.0e308]))
        @test constant_term.((a / b).coefficients) == [1, 1, 1]
        # The logarithm numerator cancels exactly before division; splitting
        # its two quotients would instead produce Inf - Inf.
        s = TimeSeries(wrap.([2.0^-1000, 2.0^10, 2.0^1019]))
        @test constant_term(coefficient(log(s), 2)) == 0
    end
    for scale in (1.0e308, 1.0e-308)
        result = TimeSeries([scale * (1 + x), scale * (1 + x), zero(x)])^0.25
        expected = scale^0.25 .* ([1, 1 / 4, -3 / 32] * [1 1 / 4 -3 / 32])
        @test all(isapprox.(hcat(coefficients.(result.coefficients)...), expected; rtol = 2.0e-14, atol = 0))
    end
    result = TimeSeries([1.0e308 + zero(x), 1.0e308 + 1.0e-100x])^0.75
    @test coefficient(coefficient(result, 1), [1]) ≈ 7.5e-178

    argument = TimeSeries([0.001one(x), 0.0006x, zero(x)])
    numerator = TimeSeries([0.001one(x), zero(x), zero(x)])
    joined_log = TimeSeries([one(x), sqrt(0.12) * x, -0.06x^2])
    a, b = TimeSeries([0.06x, 0.06x]), TimeSeries([one(x), one(x)])
    seed = TimeSeries([0.01x, 1000one(x), zero(x)])
    shifted_seed = TimeSeries([1 + 0.01x, 1000one(x), zero(x)])
    exponent = TimeSeries([zero(x), 1000one(x), zero(x)])
    base = TimeSeries([1 + 0.01x, zero(x), zero(x)])
    stored = TimeSeries([1 + 0.01x + 2x^2, 3 + 0.02x + 4x^2])
    parameter = 0.01x
    large_seed = TimeSeries([3one(x), 1.0e100one(x)])
    exact_x = TaylorPolynomial{Rational{BigInt}}(x)
    exact_seed = TimeSeries([exact_x / 100, 1000one(exact_x), zero(exact_x)])
    set_coefficient_tolerance!(0.1)
    try
        @test coefficient(coefficient(log(argument), 2), [2]) ≈ -0.18
        @test coefficient(coefficient(numerator / argument, 2), [2]) ≈ 0.36
        @test coefficient(coefficient(log(joined_log), 2), [2]) ≈ -0.12
        @test coefficient(coefficient(a * b, 1), [1]) ≈ 0.12
        # Small seed/earlier-time coefficients must survive until the entire
        # public operation is finished, including coupled trig recurrences.
        for (f, expected) in (
                (exp, [0, 10, 5000]), (sin, [0, 0, -5000]),
                (cos, [0, -10, 0]), (sinh, [0, 0, 5000]), (cosh, [0, 10, 0]),
            )
            @test [coefficient(p, [1]) for p in f(seed).coefficients] ≈ expected
        end
        @test [coefficient(p, [1]) for p in sqrt(shifted_seed).coefficients] ≈ [0, -2.5, 1875]
        @test coefficient(coefficient(shifted_seed^2, 1), [1]) ≈ 20
        @test coefficient(coefficient(shifted_seed^-2, 1), [1]) ≈ 60
        @test coefficient(coefficient(base^exponent, 2), [2]) ≈ 50
        recorded = only(taylor_expand((u, p, t) -> [exp(p + 1000t)], [zero(x)], 0.0; order = 3, parameters = parameter))
        @test [coefficient(p, [1]) for p in recorded.coefficients] ≈ [0, 0, 5, 10000 / 6]
        with_order(1) do
            @test all(p -> coefficients(p) == [constant_term(p), 0, 0], (+stored).coefficients)
            @test all(p -> coefficients(p) == [constant_term(p), 0, 0], (stored + 0).coefficients)
            @test [coefficient(p, [2]) for p in copy(stored).coefficients] == [2, 4]
        end
        set_coefficient_tolerance!(2.0)
        @test constant_term.(exp(large_seed).coefficients) ≈ exp(3) .* [1, 1.0e100]
        @test [coefficient(p, [1]) for p in exp(exact_seed).coefficients] ≈ [0, 10, 5000]
        recorded = only(taylor_expand((u, p, t) -> [3 + 10000 * (0.001t)], [zero(x)], 0.0; order = 2))
        @test constant_term.(recorded.coefficients) ≈ [0, 3, 5]
    finally
        set_coefficient_tolerance!(0.0)
    end
end

@testset "Time degree factors retain interval guarantees" begin
    IntervalArithmetic.configure(rounding = :correct)
    value = TimeSeries(interval.([2, 1, 0]))(1)
    @test IntervalArithmetic.isguaranteed(value) && IntervalArithmetic.in_interval(3, value)
    x, = variables(Interval{Float64}, 1; order = 2)
    for input in (
                TimeSeries(interval.([2, 1, 0])),
                TimeSeries([interval(2) + x, one(x), zero(x)]),
            ), f in (exp, log, sqrt, sin, t -> t^(1 // 2), integrate, differentiate)
        @test all(f(input).coefficients) do coefficient
            values = coefficient isa TaylorPolynomial ? coefficients(coefficient) : (coefficient,)
            all(c -> IntervalArithmetic.isguaranteed(c) && IntervalArithmetic.decoration(c) == IntervalArithmetic.com, values)
        end
    end
    solution = only(taylor_expand((u, p, t) -> [u[1] + 1], [interval(1)], interval(0); order = 3))
    @test all(zip(solution.coefficients, (1, 2, 1, 1 // 3))) do (actual, expected)
        IntervalArithmetic.isguaranteed(actual) && IntervalArithmetic.in_interval(expected, actual)
    end
end

@testset "Solver mutation, rejection and historical interpolation" begin
    initial, parameters = [1.0], [1.0]
    problem = ODEProblem(growth!, initial, (0.0, 1.0), parameters)
    callback = ContinuousCallback(
        (u, t, integrator) -> t - 0.5,
        integrator -> (integrator.p[1] = 2)
    )
    solution = solve(problem, TaylorMethod(12); adaptive = false, dt = 0.4, callback)
    @test SciMLBase.successful_retcode(solution)
    @test initial == [1] # No accidental alias of the caller's initial state.
    @test only(solution(0.45)) ≈ exp(0.45) rtol = 1.0e-12
    @test only(solution(0.75)) ≈ exp(1) rtol = 1.0e-12
    out = [NaN]
    solution(out, 0.75, Val{1})
    @test only(out) ≈ 2exp(1) rtol = 1.0e-12
    parameters[1] = 1
    rejected = solve(
        problem, TaylorMethod(12); dt = 1.0, abstol = 1.0e-13, reltol = 1.0e-13
    )
    @test rejected.stats.nreject > 0
    @test only(last(rejected.u)) ≈ exp(1) rtol = 1.0e-10
    integrator = init(problem, TaylorMethod(10); adaptive = false, dt = 0.1)
    step!(integrator)
    integrator.u[1] = 3
    integrator.p[1] = 2
    SciMLBase.derivative_discontinuity!(integrator, true)
    step!(integrator)
    @test only(integrator.u) ≈ 3exp(0.2) rtol = 1.0e-12
    reinit!(integrator, [3.0])
    @test integrator.t == 0
    aliased = [1.0]
    integrator = init(
        ODEProblem(growth!, aliased, (0.0, 0.2), [1.0]), TaylorMethod(10);
        adaptive = false, dt = 0.1, alias = ODEAliasSpecifier(alias_u0 = true)
    )
    step!(integrator)
    @test integrator.u === aliased && only(aliased) ≈ exp(0.1)
    x, = variables(1; order = 2)
    @test_throws ArgumentError init(
        ODEProblem(growth!, [1.0], (0.0, 1.0), [1 + x]), TaylorMethod(6)
    )
    @test_throws ArgumentError init(
        ODEProblem(growth!, [interval(1)], (0.0, 1.0), [1]), TaylorMethod(6)
    )
    @test_throws ArgumentError init(
        ODEProblem((u, p, t) -> [p[]], [0.0], (0.0, 1.0), Ref(interval(1, 2))), TaylorMethod(3)
    )
    scaled = solve(
        ODEProblem((u, p, t) -> [3.0e-320t^2], [0.0], (0.0, 1.0e110)), TaylorMethod(3);
        dt = 1.0e110, abstol = 1.0e12, reltol = 0, maxiters = 2
    )
    @test SciMLBase.successful_retcode(scaled) && scaled.stats.nreject == 0
    @test only(last(scaled.u)) ≈ Float64(BigFloat(3.0e-320) * big"1e330" / 3) rtol = 1.0e-14
end

@testset "Dense derivatives from an independent polynomial" begin
    forcing(u, p, t) = [2 + 6t + 12t^2, -4t]
    solution = solve(
        ODEProblem(forcing, [1.0, 5.0], (0.0, 0.6)), TaylorMethod(3);
        adaptive = false, dt = 0.6
    )
    expected = ([1.978, 4.82], [4.88, -1.2], [13.2, -4], [24.0, 0.0], [0.0, 0.0])
    for d in 0:4
        out = fill(NaN, 2)
        solution(out, 0.3, Val{d}; idxs = [2, 1])
        @test out ≈ reverse(expected[d + 1])
        @test solution(0.3, Val{d}; idxs = 2) ≈ expected[d + 1][2]
    end
    backward = solve(
        ODEProblem(forcing, [4.144, 4.28], (0.6, 0.0)), TaylorMethod(3);
        adaptive = false, dt = -0.6, saveat = [0.45, 0.15],
        save_start = false, save_end = false
    )
    @test backward.t == [0.45, 0.15]
    @test hcat(backward.u...) ≈ [2.872 1.381; 4.595 4.955]
    for iip in (false, true)
        forcing!(du, u, p, t) = (du .= forcing(u, p, t))
        problem = ODEProblem(iip ? forcing! : forcing, [1.0, 5.0], (0.0, 0.6))
        selected = solve(problem, TaylorMethod(3); adaptive = false, dt = 0.3, save_idxs = 2)
        @test selected(0.15) isa Real && selected(0.15) ≈ 4.955
        @test selected(0.15, Val{1}) ≈ -0.6
        @test selected([0.15, 0.45]).u ≈ [4.955, 4.595]
        integrator = init(problem, TaylorMethod(3); adaptive = false, dt = 0.3, calck = false)
        @test SciMLBase.get_du(integrator) ≈ [2, 0]
        step!(integrator)
        derivative = zeros(2)
        SciMLBase.get_du!(derivative, integrator)
        @test derivative ≈ expected[2]
        @test integrator(0.15, Val{2}) ≈ [9.6, -4]
    end
end

@testset "Flow callbacks retain independent returned states" begin
    scratch = Real[0.0]
    initial(p) = (scratch[1] = p[1]; scratch)
    checkpoint(u, span) = (scratch[1] = 2u[1]; scratch)
    for (advance, initialize, expected) in (
            ((u, span) -> u, initial, 0.25),
            (checkpoint, p -> [p[1]], 0.5),
        )
        flow = adaptive_flow(advance, initialize, [0.0], [1.0], [0.0, 1.0]; order = 2, atol = 1.0e-10, max_patches = 1)
        @test flow.converged && only(flow([0.25])) ≈ expected
    end
    function online(u, span, monitor)
        value = [u[1] + last(span) - first(span)]
        monitor(value, last(span))
        scratch[1] = value[1]
        return (; state = scratch, time = last(span))
    end
    flow = adaptive_flow(online, initial, [0.0], [1.0], (0.0, 1.0); order = 2, atol = 1.0e-10, max_patches = 1)
    @test flow.converged && only(flow([0.25])) ≈ 1.25
end
