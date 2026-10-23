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
end
