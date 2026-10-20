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
    println(
        (;
            julia = VERSION, revision, cpu = Sys.CPU_NAME,
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
    compilation = hasproperty(first_call, :compile_time) ? first_call.compile_time : missing
    return (; warmed..., first_seconds = first_call.time, compilation)
end
