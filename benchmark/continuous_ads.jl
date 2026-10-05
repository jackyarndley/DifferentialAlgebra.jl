# Optional benchmark, excluded from the test suite.
# julia --project=examples benchmark/continuous_ads.jl
using DifferentialAlgebra, IntervalArithmetic, ForwardDiff

function measure_continuity(f; samples = 10)
    f()
    times = [@elapsed f() for _ in 1:samples]
    bytes = @allocated f()
    return (; seconds = minimum(times), bytes)
end

for (label, f, dimensions, splitter) in (
        ("Diagonal nonlinear map", v -> exp(v[1] + v[2]) + (v[1] - v[2])^2, 2, :tail),
        ("Diagonal nonlinear map", v -> exp(v[1] + v[2]) + (v[1] - v[2])^2, 2, :oriented),
        ("Six-variable uncertainty map", v -> exp(sum(v) / 6) * cos(v[1] * v[2]) + sum(x -> x^2, v), 6, :tail),
    )
    for budget in (0, 32 * 1024^2)
        # Six-dimensional maps use a representative small uncertainty box.
        # A broad box can exhaust the patch budget with this simple bounder.
        radius = dimensions == 6 ? 1 // 100 : 1 // 4
        box = fill(interval(Float64, -radius, radius), dimensions)
        build() = adaptive_map(f, box; estimator = IntervalBound(), splitter, order = 3, atol = 1 // 100000, table_bytes = budget)
        source = build()
        blend() = continuous_map(f, source; continuity = :c2, order = 4, table_bytes = budget)
        s = blend()
        point = fill(Float64(radius) / 4, dimensions)
        println((; label, radius, budget, splitter, leaves = length(source.patches), source_construction = measure_continuity(build), overlap_construction = measure_continuity(blend), source_query = measure_continuity(() -> source(point)), smooth_query = measure_continuity(() -> s(point)), gradient = measure_continuity(() -> ForwardDiff.gradient(s, point)), hessian = measure_continuity(() -> ForwardDiff.hessian(s, point)), original_enclosure = measure_continuity(() -> enclose(s)), full_width = diam(enclose(s)), point_width = diam(enclose(s, point)), uniform_blend_error = sup(abs(only(s.error_bounds)))))
    end
end
