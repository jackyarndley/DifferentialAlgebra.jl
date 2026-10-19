# Optional benchmark, excluded from the test suite.
# julia --project=benchmark benchmark/continuous_ads.jl
using DifferentialAlgebra, IntervalArithmetic, ForwardDiff

include("measurement.jl")

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
        println((; label, radius, budget, splitter, patches = length(source.patches), source_construction = measure(build), overlap_construction = measure(blend), source_query = measure(() -> source(point)), smooth_query = measure(() -> s(point)), gradient = measure(() -> ForwardDiff.gradient(s, point)), hessian = measure(() -> ForwardDiff.hessian(s, point)), original_enclosure = measure(() -> enclose(s)), full_width = diam(enclose(s)), point_width = diam(enclose(s, point)), uniform_blend_error = sup(abs(only(s.error_bounds)))))
    end
end
