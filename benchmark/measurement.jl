# Shared setup for standalone scripts. Compilation is reported separately from
# warmed batched minima; run paired revisions sequentially on the same machine.
function measure(f::F; samples = 30, batch = 10) where {F}
    first_call = @timed f()
    f()
    seconds = minimum(1:samples) do _
        elapsed = @elapsed for _ in 1:batch
            f()
        end
        elapsed / batch
    end
    compilation = hasproperty(first_call, :compile_time) ? first_call.compile_time : missing
    return (; seconds, bytes = @allocated(f()), first_seconds = first_call.time, compilation)
end
