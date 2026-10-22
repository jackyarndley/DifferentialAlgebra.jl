using Test, DifferentialAlgebra, Distributed, Logging
const DA = DifferentialAlgebra

@testset "Owning contexts, inherited policies and exception safety" begin
    x, = variables(1; order = 4)
    p = x^4 + 1.0e-4x
    set_coefficient_tolerance!(0.001)
    with_order(2) do
        @test fetch(Threads.@spawn degree(p + 0)) == 0
        @test fetch(@async truncation_order()) == 2
    end
    @test_throws ErrorException with_order(1) do
        error("restore working order")
    end
    @test truncation_order() == 4
    set_coefficient_tolerance!(0)
    DA.with_algebra(3, 2) do _
        u = variable(1)
        Logging.with_logger(NullLogger()) do
            @test fetch(Threads.@spawn variable(1).algebra === u.algebra)
            @test fetch(
                Threads.@spawn (u + one(TaylorPolynomial{Float64})).algebra === u.algebra
            )
        end
        @test nvariables(p) == 1 && degree(p + 0) == 4
    end
    @test variable(1).algebra === x.algebra
    # Independent builders must neither replace nor lock the caller's algebra.
    tasks = [
        Threads.@spawn DA.with_algebra(2, n) do _
            yield()
            nvariables() == n && degree(variable(1)^2) == 2
        end for n in (2, 3)
    ]
    @test all(fetch, tasks)
    @test nvariables() == 1
    entered, release = Channel{Nothing}(1), Channel{Nothing}(1)
    builder = Threads.@spawn DA.with_algebra(2, 2) do _
        u = variable(1)
        put!(entered, nothing)
        take!(release)
        degree(u^2) == 2 && variable(1).algebra === u.algebra
    end
    bind(entered, builder)
    take!(entered)
    initialize!(3, 3)
    put!(release, nothing)
    @test fetch(builder)
    @test nvariables() == 3
end

@testset "Distributed policy persists across remote-call tasks" begin
    inherited = Base.julia_cmd().exec[2:end]
    flags = `$inherited --startup-file=no --project=$(Base.active_project())`
    workers = addprocs(1; exeflags = flags)
    worker = only(workers)
    try
        remotecall_fetch(
            Core.eval, worker, Main, quote
                using DifferentialAlgebra
                variables(1; order = 4)
                set_truncation_order!(2)
                set_coefficient_tolerance!(0.001)
            end
        )
        @test remotecall_fetch(truncation_order, worker) == 2
        @test remotecall_fetch(coefficient_tolerance, worker) == 0.001
    finally
        rmprocs(workers)
    end
end
