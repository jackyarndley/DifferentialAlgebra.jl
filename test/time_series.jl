@testset "Time series arithmetic" begin
    # Compare independent univariate and multivariate coefficient recurrences.
    x, = variables(1; order = 8)
    a = TimeSeries([1.2, 0.3, -0.1, zeros(6)...])
    b = TimeSeries([0.8, -0.2, zeros(7)...])
    ap, bp = 1.2 + 0.3x - 0.1x^2, 0.8 - 0.2x
    functions = (
        identity, -, exp, log, sqrt, cbrt, sin, cos, tan, sinh, cosh, tanh, inv,
        s -> s^3, s -> s^0, s -> s^(-3), s -> s^(1 // 3), s -> s^1.7,
        s -> s + 2, s -> 2 + s, s -> s - 2, s -> 2 - s,
        s -> 2s, s -> s * 2, s -> s / 2, s -> 2 / s, s -> 2.0^s,
    )
    for f in functions
        result, expected = f(a), f(ap)
        @test max_order(result) == 8
        @test result.coefficients ≈ [coefficient(expected, [k]) for k in 0:8] atol = 2.0e-14
    end
    for op in (+, -, *, /, ^)
        result, expected = op(a, b), op(ap, bp)
        @test result.coefficients ≈ [coefficient(expected, [k]) for k in 0:8] atol = 2.0e-14
    end
    @test differentiate(integrate(a)).coefficients ≈ a.coefficients
    @test max_order(integrate(a)) == 9
    @test max_order(differentiate(a)) == 7
    @test a(0.2) ≈ ap(0.2)
    @test evaluate(a, 0.2) == a(0.2)
    @test coefficient_type(a) == Float64
    @test iszero(zero(a)) && isfinite(a)
    @test coefficient(one(a), 0) == 1
    @test iszero(differentiate(TimeSeries([3.0])))
    @test sincos(a)[1].coefficients == sin(a).coefficients
    @test abs2(a).coefficients == (a * a).coefficients
    @test real(a) === a && conj(a) === a
    @test !isfinite(TimeSeries([Inf]))
    @test_throws ArgumentError TimeSeries(Float64[])
    @test_throws BoundsError coefficient(a, -1)
    @test_throws BoundsError coefficient(a, 9)
    for op in (+, -, *, /, ^)
        @test_throws ArgumentError op(a, TimeSeries([1.0]))
    end
    t = TimeSeries([0.0, 1.0])
    for f in (inv, log, sqrt, cbrt, s -> s^0.5)
        @test_throws DomainError f(t)
    end
    for T in (Float32, BigFloat, Rational{BigInt})
        exact = TimeSeries(T[1, -1, 0, 0, 0])
        @test inv(exact).coefficients == ones(T, 5)
        @test coefficient_type(inv(exact)) == T
        @test coefficient_type(exact / 3) == T
    end
    # The outer constructor owns its storage, including mutable polynomials.
    values = [1 + x, 2 + x]
    s = TimeSeries(values)
    set_coefficient!(values[1], [1], 7)
    @test coefficient(coefficient(s, 0), [1]) == 1
    copied = copy(s)
    set_coefficient!(coefficient(copied, 0), [1], 3)
    @test coefficient(coefficient(s, 0), [1]) == 1
end

@testset "Taylor ODE expansions" begin
    # Nonautonomous equations, a nonzero epoch and constant derivative outputs.
    initial = [2.0, 3.0]
    solution = taylor_expand((u, p, t) -> [2t, p], initial, 1.5; order = 8, parameters = 4.0)
    @test solution[1](0.2) ≈ 2 + (1.7^2 - 1.5^2)
    @test solution[2](0.2) ≈ 3.8
    @test initial == [2, 3]
    # Returning references to input series must not affect the recurrence.
    oscillation = taylor_expand((u, p, t) -> [u[2], -u[1]], [0.0, 1.0], 0.0; order = 20)
    @test oscillation[1](0.5) ≈ sin(0.5) atol = 1.0e-15
    @test oscillation[2](-0.5) ≈ cos(0.5) atol = 1.0e-15
    linear = taylor_expand(
        (u, p, t) -> p * u, [0.0, 1.0], 0.0;
        order = 20, parameters = [0.0 1.0; -1.0 0.0]
    )
    @test [s(0.5) for s in linear] ≈ [sin(0.5), cos(0.5)] atol = 1.0e-15
    swapped = taylor_expand((u, p, t) -> reverse(u), [0.0, 1.0], 0.0; order = 12)
    @test swapped[1](0.2) ≈ sinh(0.2)
    nonlinear = only(taylor_expand((u, p, t) -> u .^ 2, [1.0], 0.0; order = 16))
    @test nonlinear.coefficients == ones(17)
    logistic = only(taylor_expand((u, p, t) -> [u[1] * (1 - u[1])], [0.3], 0.0; order = 18))
    @test logistic(0.4) ≈ 1 / (1 + (1 / 0.3 - 1) * exp(-0.4)) rtol = 1.0e-14
    # Independently check every recorded elementary recurrence against the
    # existing eager TimeSeries arithmetic, using a time-forced ODE.
    τ = TimeSeries([1.2, 1.0, zeros(10)...])
    for g in (
            exp, log, sin, cos, tan, sinh, cosh, tanh, sqrt, cbrt, inv,
            s -> s^(-3), s -> s^1.7, s -> s^(1 // 3), s -> 2.0^s,
            s -> s^s, s -> sin(s) * exp(s) / (1 + s), s -> s^2,
        )
        forcing = (u, p, t) -> [g(t)]
        expansion = only(taylor_expand(forcing, [0.0], 1.2; order = 12))
        @test expansion.coefficients ≈ integrate(g(τ)).coefficients rtol = 2.0e-13 atol = 2.0e-14
    end
    for T in (Float32, BigFloat)
        solution = taylor_expand((u, p, t) -> u, T[1], zero(T); order = 20)
        @test coefficient_type(only(solution)) == T
        @test only(solution)(T(0.1)) ≈ exp(T(0.1)) rtol = max(T(1.0e-35), 5eps(T))
    end

    δ, = variables((:δ,); order = 2)
    initial_map = (1 + δ)^2
    solution = taylor_expand((u, p, t) -> u, [initial_map], 0.0; order = 20)
    @test truncation_order() == 2 && max_order() == 2
    @test max_order(only(solution)) == 20
    @test coefficient(coefficient(only(solution), 20), [2]) ≈ 1 / factorial(big(20))
    @test coefficient_norm(only(solution)(0.1) - exp(0.1) * initial_map) < 1.0e-15
    @test coefficient_norm(initial_map - (1 + δ)^2) == 0
    # Polynomial parameters and mixed numeric/polynomial time arithmetic.
    solution = taylor_expand((u, p, t) -> [p * u[1] + t], [one(δ)], 1.0; order = 8, parameters = 1 + δ)
    @test coefficient_norm(coefficient(only(solution), 1) - (2 + δ)) == 0
    @test coefficient_norm(coefficient(only(solution), 2) - (3 + 3δ + δ^2) / 2) == 0
    @test max_order() == 2

    # Uncertain parameters must promote numeric initial states, never silently
    # convert the polynomial derivatives to their constant coefficients.
    solution = taylor_expand((u, p, t) -> [p * u[1]], [1.0], 0.0; order = 5, parameters = 1 + δ)
    @test coefficient_type(only(solution)) == TaylorPolynomial{Float64}
    @test coefficient_norm(coefficient(only(solution), 2) - (1 + δ)^2 / 2) == 0
    exact = only(taylor_expand((u, p, t) -> u, Rational{BigInt}[1], 0 // 1; order = 12))
    @test coefficient(exact, 12) == 1 // factorial(big(12))

    f(u, p, t) = u
    @test_throws ArgumentError taylor_expand(f, [1.0], 0.0; order = 0)
    @test_throws ArgumentError taylor_expand(f, Float64[], 0.0)
    @test_throws ArgumentError taylor_expand(f, [Inf], 0.0)
    @test_throws ArgumentError taylor_expand(f, [1.0], NaN)
    @test_throws ArgumentError taylor_expand(f, [1.0], δ)
    @test_throws DimensionMismatch taylor_expand((u, p, t) -> [1, 2], [1.0], 0.0)
    @test_throws DimensionMismatch taylor_expand((u, p, t) -> 1.0, [1.0], 0.0)
    @test_throws ArgumentError taylor_expand((u, p, t) -> [1im], [1.0], 0.0)
    @test_throws ArgumentError taylor_expand((u, p, t) -> [TimeSeries([1.0, 1.0])], [1.0], 0.0)
    @test_throws DomainError taylor_expand((u, p, t) -> [Inf], [1.0], 0.0)
    @test δ(0.3) == 0.3 # No algebra changes, even on failure.
end
