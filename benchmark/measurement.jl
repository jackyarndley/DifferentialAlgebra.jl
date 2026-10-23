using BenchmarkTools: BenchmarkTools, @benchmark

let
    packages = [
        string(name, "=", Base.pkgversion(getfield(@__MODULE__, name)))
            for name in (
                :DifferentialAlgebra, :BenchmarkTools, :IntervalArithmetic,
                :ForwardDiff, :SpecialFunctions, :DiffEqBase, :OrdinaryDiffEqCore,
            ) if isdefined(@__MODULE__, name)
    ]
    revision = readchomp(`git -C $(pkgdir(DifferentialAlgebra)) rev-parse HEAD`)
    working_tree_dirty = !isempty(readchomp(`git -C $(pkgdir(DifferentialAlgebra)) status --porcelain`))
    println(
        (;
            julia = VERSION, revision, working_tree_dirty, cpu = Sys.CPU_NAME,
            threads = Threads.nthreads(), packages,
        )
    )
end

function warmed_measurement(f::F, samples, batch, seconds) where {F}
    trial = if batch === nothing
        @benchmark $f() samples = samples seconds = seconds
    else
        @benchmark $f() samples = samples evals = batch seconds = seconds
    end
    best = minimum(trial)
    return (;
        seconds = best.time / 1.0e9, bytes = best.memory,
        allocations = best.allocs, evals = trial.params.evals,
    )
end

function measure(f; samples = 30, batch = nothing, seconds = 0.5)
    # Dynamic entry keeps workload compilation inside the first-call timer.
    # BenchmarkTools calibrates its warmed loop after that timer has stopped.
    first_call = @timed Base.invokelatest(f)
    warmed = Base.invokelatest(warmed_measurement, f, samples, batch, seconds)
    return (; warmed..., first_seconds = first_call.time, compilation = first_call.compile_time)
end
