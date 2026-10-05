# Optional scientific benchmark, excluded from the test suite.
# julia --project=examples benchmark/polygon_ads.jl
using DifferentialAlgebra, IntervalArithmetic

function measure_polygon(f; samples = 10)
    f()
    times = [@elapsed f() for _ in 1:samples]
    bytes = @allocated f()
    return (; seconds = minimum(times), bytes)
end

for (label, f, box, order, tolerance) in (
        ("Diagonal quadratic", v -> (v[1] + v[2])^2, fill(interval(Float64, -1, 1), 2), 1, 1 // 16),
        ("Diagonal exponential", v -> exp(v[1] + v[2]), fill(interval(Float64, -1 // 2, 1 // 2), 2), 3, 1 // 1000),
    )
    for budget in (0, 32 * 1024^2), splitter in (:tail, :oriented)
        build() = adaptive_map(f, box; estimator = IntervalBound(), splitter, order, atol = tolerance, table_bytes = budget)
        fit = build()
        println((; label, budget, splitter, leaves = length(fit.patches), construction = measure_polygon(build), enclosure = measure_polygon(() -> enclose(fit)), full_width = diam(enclose(fit)), point_width = diam(fit([1 // 4, 1 // 4])), uniform_error = maximum(p -> sup(abs(only(p.error_bounds))), fit.patches)))
    end
end
