# Shared setup for standalone scripts. Compilation is reported separately from
# warmed batched minima; run paired revisions sequentially on the same machine.
function warmed_measurement(f::F, samples, batch) where {F}
    f()
    seconds = minimum(1:samples) do _
        elapsed = @elapsed for _ in 1:batch
            f()
        end
        elapsed / batch
    end
    return (; seconds, bytes = @allocated(f()))
end

function measure(f; samples = 30, batch = 10)
    # Dynamic entry keeps workload compilation inside the first-call timer.
    # The specialized batched loop is compiled afterward, outside that timer.
    first_call = @timed Base.invokelatest(f)
    warmed = Base.invokelatest(warmed_measurement, f, samples, batch)
    compilation = hasproperty(first_call, :compile_time) ? first_call.compile_time : missing
    return (; warmed..., first_seconds = first_call.time, compilation)
end
