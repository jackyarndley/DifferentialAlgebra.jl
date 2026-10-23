# Run separately with the local package and IntervalArithmetic in the environment.
# Reports construction work and enclosure widths alongside calibrated timings.
using DifferentialAlgebra, IntervalArithmetic

include("measurement.jl")
include("support/interval_workload.jl")

function report_ads(label, build, point)
    fit = build()
    validated = fit isa PiecewiseTaylorModel
    query = () -> fit(point)
    if !validated
        out = zeros(noutputs(fit))
        normalized = zeros(nvariables(fit))
        work = zeros(degree(fit) + 1)
        query = () -> evaluate!(out, fit, point, normalized, work)
    end
    error = maximum(
        p -> validated ? maximum(sup.(abs.(p.error_bounds))) : maximum(p.error_estimate),
        fit.patches
    )
    println(
        (;
            label, patches = length(fit.patches), construction = measure(build),
            query = measure(query; samples = 100),
            snapshot_bytes = Base.summarysize(fit), error,
        )
    )
    if validated
        println(
            (;
                label, enclosure = measure(() -> enclose(fit)),
                full_width = diam.(enclose(fit)), point_width = diam.(fit(point)),
            )
        )
    end
    return fit
end

f(v) = [exp(v[1] + v[2] / 2), sin(v[1] * v[2])]
for budget in (0, 32 * 1024^2)
    println("Multiplication table budget: ", budget)
    for method in (GuardedTail(), IntervalBound())
        report_ads(
            string(method),
            () -> adaptive_map(
                f, [-0.5, -0.5], [0.5, 0.5];
                estimator = method, order = 3, atol = 1.0e-5, table_bytes = budget
            ),
            [0.125, -0.25]
        )
    end
    multi(v) = (m = uncertainty_map(v); [m + k * v[1] for k in 1:8])
    box = fill(interval(Float64, -1 // 50, 1 // 50), 6)
    report_ads(
        "six-variable ordinary",
        () -> adaptive_map(
            uncertainty_map, fill(-0.02, 6), fill(0.02, 6);
            order = 3, atol = 1.0e-6, table_bytes = budget
        ), zeros(6)
    )
    report_ads(
        "six-variable",
        () -> adaptive_map(
            uncertainty_map, box; estimator = IntervalBound(),
            order = 3, atol = 1.0e-6, table_bytes = budget
        ), zeros(6)
    )
    report_ads(
        "eight-output snapshot",
        () -> adaptive_map(
            multi, box; estimator = IntervalBound(),
            order = 3, atol = 1.0e-6, table_bytes = budget
        ), zeros(6)
    )
end

for count in (8, 32, 128, 512)
    build() = adaptive_map(
        v -> v[1]^2, [interval(-1, 1)]; estimator = IntervalBound(),
        splitter = :width, order = 1, atol = 1 // count^2
    )
    fit = build()
    println(
        (;
            scaling = count, patches = length(fit.patches),
            point = measure(() -> fit([1 // 7]); samples = 100),
            point_width = diam(fit([1 // 7])),
            error = maximum(p -> sup(abs(only(p.error_bounds))), fit.patches),
        )
    )
end
