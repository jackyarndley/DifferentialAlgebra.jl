# Run separately; this benchmark is intentionally outside the default suite.
# julia --project=examples benchmark/interval_models.jl
using DifferentialAlgebra, IntervalArithmetic, LinearAlgebra

uncertainty_map(x) = exp(x[1] + x[2] / 2) * cos(x[3] - x[4]) + log(2 + x[5] * x[6])

function measure(f; samples = 10)
    f() # Compile and warm up.
    times = [@elapsed f() for _ in 1:samples]
    bytes = @allocated f()
    return (; seconds = minimum(times), bytes)
end

function report(; order = 3, table_bytes = 32 * 1024^2)
    box = fill(interval(Float64, -1 // 50, 1 // 50), 6)
    x = taylor_models(box; order, table_bytes)
    m = uncertainty_map(x)
    println("Six-variable model construction: ", measure(() -> uncertainty_map(x)))
    println("Model enclosure: ", measure(() -> enclose(m)))
    println("Model width / remainder width: ", (diam(enclose(m)), diam(remainder(m))))
    build_ads() = validated_adaptive_map(uncertainty_map, box; order, table_bytes, atol = 1.0e-6)
    ads = build_ads()
    println("Certified six-variable box ADS: ", measure(build_ads))
    println("Certified ADS enclosure: ", measure(() -> enclose(ads)))
    println(
        "ADS patches / width / largest uniform fit error: ",
        (length(ads.patches), diam(enclose(ads)), maximum(p -> sup(abs(only(p.error_bounds))), ads.patches))
    )

    # New initialization deliberately follows all measurements of the model.
    x = variables(6; order, table_bytes)
    p = uncertainty_map(x)
    println("Ordinary polynomial construction: ", measure(() -> uncertainty_map(x)))
    println("Stored-polynomial enclosure: ", measure(() -> enclose(p, box)))
    println("Stored-polynomial width: ", diam(enclose(p, box)))
    out = zero(p)
    println("Reusable ordinary multiplication: ", measure(() -> mul!(out, p, p)))
    compiled = compile(p)
    value, work, point = zeros(1), zeros(degree(compiled) + 1), zeros(6)
    println("Reusable ordinary evaluation: ", measure(() -> evaluate!(value, compiled, point, work)))
    return nothing
end

for budget in (0, 32 * 1024^2)
    println("\nMultiplication table budget: ", budget)
    report(; table_bytes = budget)
end
