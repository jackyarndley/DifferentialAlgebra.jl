using Test, DifferentialAlgebra
using ForwardDiff

include("support/time_rhs.jl")

# Exact binomial coefficients give an independent oracle for (1 + t + x)^p.
time_binomial(p, n) =
    prod(big(p) - j for j in 0:(n - 1); init = big(1) // 1) / factorial(big(n))

@testset "Floating time-series exponent weights" begin
    for (T, precision) in ((Float32, 128), (Float64, 128), (BigFloat, 128), (BigFloat, 384))
        setprecision(precision) do
            order = 12
            a = TimeSeries(T[1, 1, zeros(T, order - 1)...])
            for p in (1 // 2, 1 // 3, -2 // 3, 1 // (typemax(Int) ÷ 2), big(1) // 3)
                value = @inferred a^p
                @test value isa TimeSeries{T}
                for k in 0:order
                    expected = T(time_binomial(p, k))
                    @test coefficient(value, k) ≈ expected rtol = 32order * eps(T)
                end
            end
            for (f, p) in ((sqrt, 1 // 2), (cbrt, 1 // 3))
                value = @inferred f(a)
                @test value isa TimeSeries{T}
                for k in 0:order
                    expected = T(time_binomial(p, k))
                    @test coefficient(value, k) ≈ expected rtol = 32order * eps(T)
                end
                @test value.coefficients !== a.coefficients
            end
            @test a.coefficients == T[1, 1, zeros(T, order - 1)...]
        end
    end
    # Integral exponents take the exact power path, including beyond Float64's
    # exact integer range. No floating weight conversion belongs on that path.
    a = TimeSeries(Rational{BigInt}[1, 1, 0, 0])
    n = big(2)^60 + 3
    @test (a^n).coefficients == [binomial(n, k) for k in 0:3]

    # A BigInt rational uses BigFloat intermediates, even with Float32 storage.
    # Narrowing its exponent first would underflow and lose this finite slope.
    p = big(1) // big(2)^150
    a = TimeSeries(Float32[1, big(2)^100, 0])
    @test coefficient(a^p, 1) == Float32(big(1) // big(2)^50)
end

@testset "Polynomial and nested dual time-series weights" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        DifferentialAlgebra.with_algebra(3, 1; table_bytes = budget) do _
            x = variable(1, T)
            a = TimeSeries([one(x) + x, one(x), [zero(x) for _ in 1:5]...])
            # The active factory differs from the stored coefficient's owner.
            DifferentialAlgebra.with_algebra(1, 2) do _
                for (f, p) in ((sqrt, 1 // 2), (cbrt, 1 // 3))
                    value = f(a)
                    @test coefficient_type(value) === TaylorPolynomial{T}
                    for k in 0:6, m in 0:3
                        expected = T(time_binomial(p, k + m) * binomial(k + m, k))
                        actual = coefficient(coefficient(value, k), [m])
                        @test actual ≈ expected rtol = 256eps(T)
                    end
                    @test coefficient(value, 0).coeffs !== coefficient(a, 0).coeffs
                end
                @test nvariables() == 2
            end
            @test nvariables() == 1
        end
    end
    c0 = ForwardDiff.Dual(
        ForwardDiff.Dual(1.0, 1.0), ForwardDiff.Dual(1.0, 0.0)
    )
    a = TimeSeries([c0, one(c0), [zero(c0) for _ in 1:5]...])
    for (f, p) in ((sqrt, 1 // 2), (cbrt, 1 // 3), (s -> s^(2 // 3), 2 // 3))
        value = f(a)
        @test coefficient_type(value) === typeof(c0)
        for k in 0:6
            c = coefficient(value, k)
            first = Float64(time_binomial(p, k + 1) * (k + 1))
            second = Float64(time_binomial(p, k + 2) * (k + 1) * (k + 2))
            @test ForwardDiff.value(ForwardDiff.value(c)) ≈ Float64(time_binomial(p, k))
            @test only(ForwardDiff.partials(ForwardDiff.value(c))) ≈ first
            @test ForwardDiff.value(only(ForwardDiff.partials(c))) ≈ first
            @test only(ForwardDiff.partials(only(ForwardDiff.partials(c)))) ≈ second
        end
    end
end

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
        @test_throws DifferentialAlgebra.TaylorError f(t)
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
    oscillation = taylor_expand(oscillator_state, [0.0, 1.0], 0.0; order = 20)
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
        solution = taylor_expand(exponential_growth, T[1], zero(T); order = 20)
        @test coefficient_type(only(solution)) == T
        @test only(solution)(T(0.1)) ≈ exp(T(0.1)) rtol = max(T(1.0e-35), 5eps(T))
    end

    δ, = variables((:δ,); order = 2)
    initial_map = (1 + δ)^2
    solution = taylor_expand(exponential_growth, [initial_map], 0.0; order = 20)
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
    solution = taylor_expand(parameter_growth, [1.0], 0.0; order = 5, parameters = 1 + δ)
    @test coefficient_type(only(solution)) == TaylorPolynomial{Float64}
    @test coefficient_norm(coefficient(only(solution), 2) - (1 + δ)^2 / 2) == 0
    exact = only(taylor_expand(exponential_growth, Rational{BigInt}[1], 0 // 1; order = 12))
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
    @test_throws DifferentialAlgebra.TaylorError taylor_expand((u, p, t) -> [Inf], [1.0], 0.0)
    @test δ(0.3) == 0.3 # No algebra changes, even on failure.
end

@testset "Recorded multiply-adds" begin
    forcing(u, p, t) = [
        muladd(t, t, t), muladd(t, t, 2), muladd(t, 2, t), muladd(2, t, t),
        muladd(t, 2, 2), muladd(2, t, 2), muladd(2, 2, t),
    ]
    # Integrate each polynomial forcing directly, independently of the time tape.
    expected = (
        [0, 0, 1 // 2, 1 // 3, 0], [0, 2, 0, 1 // 3, 0],
        [0, 0, 3 // 2, 0, 0], [0, 0, 3 // 2, 0, 0],
        [0, 2, 1, 0, 0], [0, 2, 1, 0, 0], [0, 4, 1 // 2, 0, 0],
    )
    for T in (Float64, Rational{BigInt})
        solution = taylor_expand(forcing, zeros(T, 7), zero(T); order = 4)
        @test all(s -> coefficient_type(s) == T, solution)
        @test all(s.coefficients == T.(c) for (s, c) in zip(solution, expected))
    end
end

@testset "Tangent series and recorded integrals" begin
    s = TimeSeries([0.0, 1.0, zeros(6)...])
    for (f, signs) in ((tan, ones(8)), (tanh, [(-1)^((k - 1) ÷ 2) for k in 1:8]))
        expected = [0, 1, 0, 1 / 3, 0, 2 / 15, 0, 17 / 315] .* signs
        @test f(s).coefficients ≈ expected
        expansion = only(taylor_expand((u, p, t) -> [f(t)], [0.0], 0.0; order = 8))
        @test expansion.coefficients ≈ [0; expected ./ (1:8)]
    end
end

@testset "Time-series analytic domains" begin
    for T in (Float32, Float64, BigFloat)
        t = TimeSeries(T[0, 1, 0, 0])
        for f in (log, sqrt, s -> s^(1 // 2), s -> (-2)^s)
            @test_throws DifferentialAlgebra.TaylorError f(-one(t) + t)
        end
        @test iszero(T(0)^(2 + t))
        @test coefficient(T(0)^zero(t), 0) == 1
        @test coefficient(T(-2)^(2 + zero(t)), 0) == 4
        @test_throws DifferentialAlgebra.TaylorError T(0)^t
        @test_throws DifferentialAlgebra.TaylorError T(0)^(-one(t))
        p = T(2)^t
        @test coefficient_type(p) === T
        @test coefficient(p, 1) ≈ log(T(2))
        for f in (log, sqrt, s -> s^(1 // 2))
            @test_throws DifferentialAlgebra.TaylorError taylor_expand(
                (u, p, t) -> [f(u[1])], T[-1], T(0); order = 3
            )
        end
        @test iszero(
            only(
                taylor_expand(
                    (u, p, t) -> [T(0)^u[1]], T[2], T(0); order = 3
                )
            ) - 2
        )
    end
end


@testset "Recorded derivative buffers accept scalar zero filling" begin
    for T in (Float32, Float64, BigFloat), strategy in (:fill, :broadcast)
        series = taylor_expand(
            recorded_oscillator, T[0, 1], zero(T); order = 12, parameters = strategy
        )
        step = one(T) / 10
        @test [s(step) for s in series] ≈ [sin(step), cos(step)] rtol = 16eps(T) atol =
            2step^13 / T(factorial(13))
        @test all(s -> coefficient_type(s) === T, series)
    end
end
