using Test, DifferentialAlgebra

@testset verbose = true "TEST 1: DifferentialAlgebra initialization" begin
    @testset "1.1 Maximum order, number of variables, and number of monomials" begin
        for order in 1:10
            for nvar in 1:6
                @testset "order=$(order), nvar=$(nvar)" begin
                    DifferentialAlgebra.initialize!(order, nvar)
                    @test order == DifferentialAlgebra.max_order()
                    @test nvar == DifferentialAlgebra.nvariables()
                    @test factorial(order + nvar) / (factorial(order) * factorial(nvar)) == DifferentialAlgebra.nmonomials()
                end
            end
        end
    end

    @testset "1.2 Truncation order" begin
        for k_user in 1:10
            @testset "1.2 Trunction order (k_user=$k_user)" begin
                DifferentialAlgebra.initialize!(10, 6)

                DifferentialAlgebra.set_truncation_order!(k_user)
                k_da = DifferentialAlgebra.truncation_order()
                @test k_da == k_user
            end
        end
    end

    @testset "1.3 Cutoff for the coefficients" begin
        for cutoff_user in (1.0e-14, 1.0e-15, 1.0e-16)
            @testset "1.3 Cutoff for the coefficients (cutoff_user=$cutoff_user)" begin
                DifferentialAlgebra.initialize!(10, 6)
                DifferentialAlgebra.set_coefficient_tolerance!(cutoff_user)
                cutoff_da = DifferentialAlgebra.coefficient_tolerance()
                @test cutoff_da == cutoff_user
            end
        end
    end
end


include("support/processes.jl")
@testset "Task-local temporary algebras and numerical policies" begin
    script = """
    using DifferentialAlgebra, Test, Logging, LinearAlgebra, SpecialFunctions
    using IntervalArithmetic
    const DA = DifferentialAlgebra
    @assert Threads.nthreads() == 4

    function paused(f, verify)
        entered, release = Channel{Nothing}(1), Channel{Nothing}(1)
        worker = Threads.@spawn f(() -> (put!(entered, nothing); take!(release)))
        take!(entered)
        try
            verify()
        finally
            put!(release, nothing)
            fetch(worker)
        end
    end

    x, = variables(1; order = 4)
    root = x.algebra
    quartic, small = x^4, 1e-4x
    paused(() -> begin
        @test degree(x) == 1 && nvariables() == 1
        @test variable(1).algebra === root
        @test coefficient(x^4, [4]) == 1
    end) do hold
        DA.with_algebra(3, 2) do _
            Logging.with_logger(NullLogger()) do
                @assert nvariables() == 2 && max_order() == 3
                hold()
                @assert nvariables() == 2
            end
        end
    end
    paused(() -> begin
        @test truncation_order() == 4
        @test coefficient(x^4, [4]) == 1
    end) do hold
        with_order(2) do
            hold()
            @assert truncation_order() == 2
        end
    end
    paused(() -> begin
        @test coefficient_tolerance() == 0
        @test coefficient(small + 0, [1]) == 1e-4
    end) do hold
        set_coefficient_tolerance!(0.01)
        hold()
        @assert coefficient_tolerance() == 0.01
    end

    # New tasks start with the global basis's full order and zero tolerance.
    set_truncation_order!(2)
    set_coefficient_tolerance!(0.001)
    @test fetch(Threads.@spawn begin
        @assert truncation_order() == 4 && coefficient_tolerance() == 0
        @assert coefficient(quartic + 0, [4]) == 1
        @assert coefficient(small + 0, [1]) == 1e-4
        true
    end)
    @test truncation_order() == 2 && coefficient_tolerance() == 0.001
    @test_throws ErrorException with_order(1) do
        set_coefficient_tolerance!(0.002)
        error("restore only the working order")
    end
    @test truncation_order() == 2 && coefficient_tolerance() == 0.002
    DA.with_algebra(3, 2) do _
        # The parent's policy remains associated with its owning context.
        @test coefficient(quartic + 0, [4]) == 0
        @test coefficient(small + 0, [1]) == 0
        @test zero(x).algebra === root && one(x).algebra === root
    end
    set_truncation_order!(4)
    set_coefficient_tolerance!(0)

    setprecision(512) do
        tiny = BigFloat("1e-450") * x
        tolerance = BigFloat("1e-400")
        set_coefficient_tolerance!(tolerance)
        @test coefficient(tiny + 0, [1]) == 0
        @test fetch(Threads.@spawn coefficient(tiny + 0, [1])) ==
            coefficient(tiny, [1])
        @test coefficient_tolerance() == tolerance
    end
    set_coefficient_tolerance!(0)

    jobs = [Threads.@spawn begin
        order = mod(index, 4) + 1
        tolerance = iseven(index) ? 0.001 : 0.0
        with_order(order) do
            set_coefficient_tolerance!(tolerance)
            for _ in 1:100
                q = quartic + small
                @assert coefficient(q, [4]) == (order == 4 ? 1 : 0)
                @assert coefficient(q, [1]) == (tolerance == 0 ? 1e-4 : 0)
                yield()
            end
        end
        true
    end for index in 1:8]
    @test all(fetch, jobs)
    @test truncation_order() == 4 && coefficient_tolerance() == 0

    DA.with_algebra(4, 2) do _
        u, v = variable(1), variable(2)
        c = 2one(u)
        mgf = exp(2u + v + u^2 / 2 + 2v^2)
        @test fetch(Threads.@spawn begin
            # These objects belong to the temporary two-variable context,
            # while this task's factories still use the one-variable default.
            @assert nvariables() == 1
            constants = (zero(u), one(u), sin(c), cos(c), sinh(c),
                SpecialFunctions.besselk(1, c), zero(u)^0.5)
            @assert all(p -> p.algebra === u.algebra, constants)
            @assert collect(constant_term.(constants)) ≈
                [0, 1, sin(2), cos(2), sinh(2), SpecialFunctions.besselk(1, 2), 0]
            _, raw = DA.raw_moments(mgf, 2)
            _, central = DA.central_moments(mgf, 2)
            @assert raw ≈ [1, 2, 1, 5, 2, 5]
            @assert central ≈ [1, 0, 0, 1, 0, 4]
            values, vectors = DA.eigenpairs([1 + u zero(u); zero(u) 3 + v])
            @assert all(p -> p.algebra === u.algebra, values)
            @assert all(p -> p.algebra === u.algebra, vectors)
            true
        end)
    end

    escaped = Ref{Any}()
    @test_throws ErrorException DA.with_algebra(2, 2) do _
        escaped[] = variable(1)
        initialize!(3, 3)
        @test_throws ArgumentError degree(escaped[])
        escaped[] = variable(2)
        @test nvariables() == 3
        error("restore the caller after local reinitialization")
    end
    @test_throws ArgumentError degree(escaped[])
    @test variable(1).algebra === root && degree(x) == 1
    @test nvariables() == 1 && truncation_order() == 4

    for kind in (:box, :polygon, :interval, :continuous)
        result = Ref{Any}()
        paused(() -> begin
            for _ in 1:100
                @test variable(1).algebra === root
                @test coefficient(x^4, [4]) == 1
            end
            # A paused builder must not hold an initialization lock needed by
            # an independent task's own temporary algebra.
            other = adaptive_map(sum, zeros(3), ones(3); order = 1)
            @test other([0.1, 0.2, 0.3]) ≈ 0.6
        end) do hold
            first_call = Ref(true)
            callback(v) = begin
                if first_call[]
                    first_call[] = false
                    hold()
                end
                v[1]^2 + v[2]
            end
            a = if kind == :polygon
                polygon = ConvexPolygon([(-1,-1),(1,-1),(1,1),(-1,1)])
                adaptive_map(callback, polygon; order = 2)
            elseif kind == :interval
                adaptive_map(callback, [interval(-1,1), interval(-1,1)];
                    order = 2, estimator = IntervalBound())
            else
                adaptive_map(callback, [-1.,-1.], [1.,1.]; order = 2)
            end
            result[] = kind == :continuous ? continuous_map(callback, a) : a
        end
        value = result[]([0.25, 0.5])
        @test value isa Interval ? in_interval(9//16, value) : value ≈ 9//16
        @test nvariables() == 1 && variable(1).algebra === root
    end
    """
    @test fresh_julia(script; threads = 4, timeout = 240)
end
