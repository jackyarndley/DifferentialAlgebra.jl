using OrdinaryDiffEqCore
using SciMLBase
using SciMLBase: successful_retcode

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
