# julia --startup-file=no --project=benchmark benchmark/interval_cleanup.jl
using DifferentialAlgebra, IntervalArithmetic

include("measurement.jl")
include("support/interval_workload.jl")

function report(::Type{T}; table_bytes = 0) where {T}
    println("\nEndpoint type / table budget: ", (T, table_bytes))
    box = fill(interval(T, -1 // 50, 1 // 50), 6)
    x = taylor_models(box; order = 3, table_bytes)
    m = uncertainty_map(x)
    snapshot = compile(m)
    subbox = fill(interval(T, -1 // 100, 1 // 100), 6)
    println("Model construction: ", measure(() -> uncertainty_map(x)))
    for (label, operation) in (
            ("Model + scalar", () -> m + 1 // 3),
            ("Scalar + model", () -> 1 // 3 + m),
            ("Model - scalar", () -> m - 1 // 3),
            ("Scalar - model", () -> 1 // 3 - m),
        )
        result = operation()
        println(label, ": ", measure(operation; samples = 1000, seconds = 1.0))
        println(
            label, " range / remainder: ",
            (bound_summary(enclose(result)), bound_summary(remainder(result)))
        )
    end
    for (label, query) in (
            ("Live enclosure", () -> enclose(m)),
            ("Snapshot enclosure", () -> enclose(snapshot)),
            ("Live subbox", () -> enclose(m, subbox)),
            ("Snapshot subbox", () -> enclose(snapshot, subbox)),
        )
        println(label, ": ", measure(query; samples = 1000, seconds = 1.0))
        println(label, " bounds: ", bound_summary(query()))
    end
    return nothing
end

for budget in (0, 32 * 1024^2)
    report(Float64; table_bytes = budget)
    "--bigfloat" in ARGS && report(BigFloat; table_bytes = budget)
end
