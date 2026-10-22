using Test, DifferentialAlgebra
using OrdinaryDiffEqCore, SciMLBase
using SciMLBase: successful_retcode

import DiffEqBase
using IntervalArithmetic: Interval, interval

include("support/time_rhs.jl")

@testset "In-place zero filling and parameter validation" begin
    for strategy in (:fill, :broadcast)
        initial = [0.0, 1.0]
        solution = solve(
            ODEProblem(zeroed_oscillator!, initial, (0.0, 1.0), strategy),
            TaylorMethod(12); adaptive = false, dt = 0.1
        )
        @test successful_retcode(solution)
        @test last(solution.u) ≈ [sin(1.0), cos(1.0)] rtol = 1.0e-13
        @test initial == [0, 1]
    end
    x, = variables(1; order = 2)
    for rhs in (parameter_growth, parameter_growth!)
        for p in (1 + x, [1 + x], (; rate = 1 + x))
            problem = ODEProblem(rhs, [1.0], (0.0, 1.0), p)
            @test_throws r"Use polynomial initial states" init(problem, TaylorMethod(12))
        end
        for p in (interval(1), [interval(1)], (; rate = interval(1)))
            problem = ODEProblem(rhs, [1.0], (0.0, 1.0), p)
            @test_throws r"interval parameters" init(problem, TaylorMethod(12))
        end
    end
end

@testset "SciML Taylor integration" begin
    oscillator(u, p, t) = [u[2], -u[1]]
    function oscillator!(du, u, p, t)
        du[1], du[2] = u[2], -u[1]
        return nothing
    end
    @test_throws ArgumentError TaylorMethod(2)
    @test string(TaylorMethod(18)) == "TaylorMethod(18)"
    for f in (oscillator, oscillator!)
        problem = ODEProblem(f, [0.0, 1.0], (0.0, 2π))
        solution = solve(problem, TaylorMethod(18); abstol = 1.0e-12, reltol = 1.0e-12)
        @test successful_retcode(solution)
        @test last(solution.u) ≈ [0, 1] atol = 2.0e-11
        for t in (0.1, 1.1, 2.8, 5.7)
            @test solution(t) ≈ [sin(t), cos(t)] atol = 2.0e-11
            @test solution(t, Val{1}) ≈ [cos(t), -sin(t)] atol = 3.0e-11
            @test solution(t; idxs = 1) ≈ sin(t) atol = 2.0e-11
            @test solution(t; idxs = [2, 1]) ≈ [cos(t), sin(t)] atol = 2.0e-11
        end
        sampled = solve(problem, TaylorMethod(18); abstol = 1.0e-12, reltol = 1.0e-12, saveat = 0.1, dense = false)
        @test maximum(maximum(abs, u - [sin(t), cos(t)]) for (t, u) in zip(sampled.t, sampled.u)) < 2.0e-11
        backward = solve(remake(problem; tspan = (0.0, -2π)), TaylorMethod(18); abstol = 1.0e-12, reltol = 1.0e-12)
        @test last(backward.u) ≈ [0, 1] atol = 2.0e-11
        @test backward(-0.3) ≈ [sin(-0.3), cos(-0.3)] atol = 2.0e-11
        fixed = solve(problem, TaylorMethod(12); adaptive = false, dt = 0.1, save_everystep = false)
        @test last(fixed.u) ≈ [0, 1] atol = 2.0e-11
        iterator = init(problem, TaylorMethod(14); abstol = 1.0e-12, reltol = 1.0e-12)
        step!(iterator, 0.5, true)
        @test iterator.u ≈ [sin(0.5), cos(0.5)] atol = 2.0e-11
        @test get_du(iterator) ≈ [cos(0.5), -sin(0.5)] atol = 2.0e-11
    end

    # Constant derivatives in an in-place RHS and array tolerances.
    function forced!(du, u, p, t)
        du[1] = 2t
        du[2] = 3.0
        return nothing
    end
    solution = solve(
        ODEProblem(forced!, [1.0, 2.0], (1.0, 2.0)), TaylorMethod(8);
        abstol = [1.0e-12, 2.0e-12], reltol = [1.0e-12, 2.0e-12], dtmax = 0.2
    )
    @test last(solution.u) ≈ [4, 5]
    @test maximum(diff(solution.t)) <= 0.2 + 2eps()

    # Callback changes must be reflected in the next expansion; earlier dense
    # output must retain the expansion of its own step.
    condition(u, t, integrator) = t - 0.5
    affect!(integrator) = (integrator.u[1] *= 2)
    callback = ContinuousCallback(condition, affect!; save_positions = (true, true))
    problem = ODEProblem((u, p, t) -> u, [1.0], (0.0, 1.0))
    solution = solve(problem, TaylorMethod(16); abstol = 1.0e-12, reltol = 1.0e-12, callback)
    @test last(solution.u)[1] ≈ 2exp(1.0) rtol = 2.0e-11
    @test solution(0.2)[1] ≈ exp(0.2) rtol = 2.0e-11
    @test solution(0.8)[1] ≈ 2exp(0.8) rtol = 2.0e-11
    stop = DiscreteCallback((u, t, integrator) -> t >= 0.5, terminate!)
    stopped = solve(problem, TaylorMethod(16); callback = stop, tstops = [0.5])
    @test stopped.retcode == ReturnCode.Terminated
    @test last(stopped.t) == 0.5

    # Parameter changes can also change the recorded graph's size.
    switched(u, p, t) = iszero(p[1]) ? [zero(u[1])] : [p[1] * u[1] + 0sin(t)]
    switch = DiscreteCallback((u, t, integrator) -> t == 0.5, integrator -> (integrator.p[1] = 0.0))
    changed = solve(
        ODEProblem(switched, [1.0], (0.0, 1.0), [1.0]), TaylorMethod(16);
        callback = switch, tstops = [0.5], abstol = 1.0e-12, reltol = 1.0e-12
    )
    @test last(changed.u)[1] ≈ exp(0.5) rtol = 2.0e-11
    @test changed(0.2)[1] ≈ exp(0.2) rtol = 2.0e-11
    high = solve(problem, TaylorMethod(28); adaptive = false, dt = 0.25)
    @test high(0.0, Val{25})[1] ≈ 1 rtol = 1.0e-13
    for T in (Float32, BigFloat)
        tolerance = T === Float32 ? T(1.0e-6) : T(1.0e-40)
        precise = solve(
            ODEProblem((u, p, t) -> u, T[1], (zero(T), one(T))), TaylorMethod(20);
            abstol = tolerance, reltol = tolerance
        )
        @test last(precise.u) isa Vector{T}
        @test last(precise.u)[1] ≈ exp(one(T)) rtol = 10tolerance
    end

    δ, = variables((:δ,); order = 2)
    initial = [1 + δ + δ^2]
    map_problem = ODEProblem((u, p, t) -> u, initial, (0.0, 1.0))
    map_solution = solve(map_problem, TaylorMethod(16); abstol = 1.0e-12, reltol = 1.0e-12)
    @test coefficient_norm(last(map_solution.u)[1] - exp(1.0) * initial[1]) < 2.0e-11
    @test coefficient_norm(map_solution(0.3)[1] - exp(0.3) * initial[1]) < 2.0e-11
    @test coefficient_norm(initial[1] - (1 + δ + δ^2)) == 0
    @test max_order() == 2
    @test_throws ArgumentError solve(ODEProblem((u, p, t) -> [p * u[1]], [1.0], (0.0, 1.0), 1 + δ), TaylorMethod(12))
end

@testset "Public dense interpolation derivatives" begin
    # Independent polynomials: 1+2t+3t²+4t³ and 5-2t².
    rhs(u, p, t) = [2 + 6t + 12t^2, -4t]
    problem = ODEProblem(rhs, [1.0, 5.0], (0.0, 0.6))
    solution = solve(problem, TaylorMethod(3); adaptive = false, dt = 0.6)
    expected = ([1.978, 4.82], [4.88, -1.2], [13.2, -4], [24.0, 0.0], [0.0, 0.0])
    for D in 0:4, idxs in (nothing, [2, 1], 2, view([1, 2], 1:1))
        indices = idxs === nothing ? [1, 2] : idxs isa Integer ? [idxs] : idxs
        out = fill(NaN, length(indices))
        wanted = expected[D + 1][indices]
        solution(out, 0.3, Val{D}; idxs)
        @test out ≈ wanted
        actual = solution(0.3, Val{D}; idxs)
        @test (idxs isa Integer ? [actual] : actual) ≈ wanted
    end
    @test_throws DimensionMismatch solution(zeros(1), 0.3)
    @test_throws ArgumentError solution(zeros(2), 0.3, Val{-1})
end

@testset "TaylorMethod rejects interval coefficient states clearly" begin
    x, = variables(Interval{Float64}, 1; order = 2)
    for initial in ([1 + x], [interval(1)])
        problem = ODEProblem((u, p, t) -> u, initial, (0.0, 1.0))
        @test_throws ArgumentError init(problem, TaylorMethod(8))
    end
end

@testset "In-place Taylor state ownership and initialized derivatives" begin
    concrete = Ref(true)
    function growth!(du, u, p, t)
        concrete[] &= isconcretetype(eltype(du))
        du[1] += u[1]
        # The untouched second derivative remains zero.
        return nothing
    end
    for alias in (false, true)
        initial = [1.0, 2.0]
        problem = ODEProblem(growth!, initial, (0.0, 0.3))
        integrator = init(
            problem, TaylorMethod(12); adaptive = false, dt = 0.1,
            alias = ODEAliasSpecifier(alias_u0 = alias)
        )
        storage = integrator.u
        @test (storage === initial) == alias
        @test integrator.uprev !== storage
        for _ in 1:3
            step!(integrator)
            @test integrator.u === storage
            @test integrator.u ≈ [exp(integrator.t), 2] rtol = 1.0e-13
        end
        @test concrete[]
        @test initial ≈ (alias ? [exp(0.3), 2] : [1, 2]) rtol = 1.0e-13
        @test integrator.sol(0.05) ≈ [exp(0.05), 2] rtol = 1.0e-13
    end
end

@testset "In-place callbacks preserve aliased state storage" begin
    growth!(du, u, p, t) = (du[1] = u[1]; nothing)
    initial = [1.0]
    callback = ContinuousCallback(
        (u, t, integrator) -> t - 0.15,
        integrator -> (integrator.u[1] *= 2)
    )
    integrator = init(
        ODEProblem(growth!, initial, (0.0, 0.3)), TaylorMethod(14);
        adaptive = false, dt = 0.1, callback,
        alias = ODEAliasSpecifier(alias_u0 = true)
    )
    solution = solve!(integrator)
    @test integrator.u === initial
    @test only(initial) ≈ 2exp(0.3) rtol = 1.0e-13
    @test only(solution(0.05)) ≈ exp(0.05) rtol = 1.0e-13
    @test only(solution(0.25)) ≈ 2exp(0.25) rtol = 1.0e-13
end

@testset "In-place coefficient types and polynomial ownership" begin
    growth!(du, u, p, t) = (du[1] = u[1]; nothing)
    for T in (Float32, BigFloat)
        initial = T[1]
        endpoint = T(3) / 10
        solution = solve(
            ODEProblem(growth!, initial, (zero(T), endpoint)),
            TaylorMethod(14); adaptive = false, dt = one(T) / 10
        )
        @test eltype(last(solution.u)) === T
        @test only(last(solution.u)) ≈ exp(endpoint) rtol =
            (T === Float32 ? T(1.0e-6) : T(1.0e-24))
        @test initial == T[1]
    end
    x, = variables(1; order = 2)
    polynomial = 1 + x
    initial = [polynomial]
    integrator = init(
        ODEProblem(growth!, initial, (0.0, 0.3)), TaylorMethod(14);
        adaptive = false, dt = 0.1, alias = ODEAliasSpecifier(alias_u0 = true)
    )
    solve!(integrator)
    @test integrator.u === initial
    @test coefficient_norm(only(initial) - exp(0.3) * (1 + x)) < 1.0e-13
    @test coefficient_norm(polynomial - (1 + x)) == 0
end

@testset "Rejected steps reuse the expansion" begin
    for inplace in (false, true)
        calls = Ref(0)
        growth(u, p, t) = (calls[] += 1; u)
        growth!(du, u, p, t) = (calls[] += 1; du[1] = u[1]; nothing)
        problem = ODEProblem(inplace ? growth! : growth, [1.0], (0.0, 1.0))
        solution = solve(
            problem, TaylorMethod(8); dt = 1.0,
            abstol = 1.0e-14, reltol = 1.0e-14
        )
        @test solution.stats.nreject > 0
        # Numeric promotion and the initial graph take two calls. Each accepted
        # step then needs one graph; rejected dt retries do not call the RHS.
        @test calls[] == solution.stats.naccept + 2
        @test solution.stats.nf == calls[]
        @test only(last(solution.u)) ≈ exp(1.0) rtol = 1.0e-12
    end
end

@testset "Expansion invalidation and historical parameter snapshots" begin
    for inplace in (false, true)
        growth(u, p, t) = p[1] .* u
        growth!(du, u, p, t) = (du[1] = p[1] * u[1]; nothing)
        rhs = inplace ? growth! : growth
        callback = ContinuousCallback(
            (u, t, integrator) -> t - 0.5,
            integrator -> (integrator.p[1] = 2.0)
        )
        solution = solve(
            ODEProblem(rhs, [1.0], (0.0, 1.0), [1.0]), TaylorMethod(18);
            adaptive = false, dt = 0.4, callback
        )
        @test only(solution(0.45)) ≈ exp(0.45) rtol = 1.0e-13
        @test only(solution(0.75)) ≈ exp(1.0) rtol = 1.0e-13
        @test only(solution(0.45, Val{1})) ≈ exp(0.45) rtol = 1.0e-13
        @test only(solution(0.75, Val{1})) ≈ 2exp(1.0) rtol = 1.0e-13
        for reinit_cache in (false, true)
            integrator = init(
                ODEProblem(rhs, [1.0], (0.0, 0.2), [1.0]),
                TaylorMethod(14); adaptive = false, dt = 0.1
            )
            step!(integrator)
            integrator.u[1] = 3
            integrator.p[1] = 2
            SciMLBase.derivative_discontinuity!(integrator, true)
            step!(integrator)
            @test only(integrator.u) ≈ 3exp(0.2) rtol = 1.0e-13
            reinit!(integrator, [3.0]; reinit_cache)
            solution = solve!(integrator)
            @test only(last(solution.u)) ≈ 3exp(0.4) rtol = 1.0e-13
        end
    end
end
