# julia --startup-file=no --project=benchmark benchmark/interval_kernels.jl --profile
using DifferentialAlgebra, IntervalArithmetic
using Profile

include("measurement.jl")
include("support/interval_workload.jl")

const DA = DifferentialAlgebra
const IA = IntervalArithmetic
const EXT = Base.get_extension(DA, :DifferentialAlgebraIntervalArithmeticExt)

# Longer trials reduce clock and scheduling variation on heterogeneous CPUs.
measure_kernel(f) = measure(f; samples = 1000, seconds = 1.0)

function dense_models(; order = 3, table_bytes = 0)
    x = taylor_models(fill(interval(-1, 1), 6); order, table_bytes)
    ctx = first(x)._polynomial.algebra
    count = ctx.basis.ends[order + 1]
    p, q = DA.allocate(ctx, Interval{Float64}, count),
        DA.allocate(ctx, Interval{Float64}, count)
    for k in eachindex(p.coeffs)
        p.coeffs[k] = interval(Float64, (-1)^k * (k % 5 + 1) // 1000)
        q.coeffs[k] = interval(Float64, (k % 7 + 1) // 2000)
    end
    p.len = q.len = length(p.coeffs)
    return TaylorModel(p, interval(Float64, -1 // 10000, 1 // 10000), first(x)),
        TaylorModel(q, interval(Float64, -1 // 10000, 1 // 10000), first(x))
end

function reference_power_cache(box, order, ::Type{I}) where {I}
    powers = Matrix{I}(undef, length(box), order + 1)
    for v in eachindex(box), k in 0:order
        powers[v, k + 1] = IA.pown(box[v], k)
    end
    return powers
end

function profile_workload(f::F; seconds = 2) where {F}
    f()
    Profile.clear()
    deadline = time_ns() + seconds * 10^9
    Profile.@profile while time_ns() < deadline
        for _ in 1:50
            f()
        end
    end
    Profile.print(; format = :flat, sortedby = :count, mincount = 20)
    return nothing
end

function report(; table_bytes = 0, profiling = false)
    println("\nMultiplication table budget: ", table_bytes)
    a, b = dense_models(; table_bytes)
    p, q = a._polynomial, b._polynomial
    # Measure first-call product compilation before compiling its subkernels.
    println("Dense model product: ", measure_kernel(() -> a * b))
    powers = EXT.power_cache(a._coordinates.normalized, 6, Interval{Float64})
    bounds = EXT.unit_product_bounds(powers, 6, Interval{Float64})
    println(
        "Normalized powers, reference: ", measure_kernel(
            () ->
            reference_power_cache(a._coordinates.normalized, 6, Interval{Float64})
        )
    )
    println(
        "Normalized powers, kernel: ", measure_kernel(
            () ->
            EXT.power_cache(a._coordinates.normalized, 6, Interval{Float64})
        )
    )
    if isdefined(EXT, :decorated_discarded_unit_bound)
        println(
            "Discarded products, decorated reference: ", measure_kernel(
                () ->
                EXT.decorated_discarded_unit_bound(p, q, 3, bounds, Interval{Float64})
            )
        )
    end
    println(
        "Dense polynomial range: ", measure_kernel(
            () -> EXT.polynomial_bound(
                p, powers,
                Interval{Float64}
            )
        )
    )
    println(
        "Discarded products: ", measure_kernel(
            () -> EXT.discarded_bound(
                p, q, 3,
                powers, Interval{Float64}
            )
        )
    )
    println("Dense model range: ", bound_summary(enclose(a * b)))
    println("Dense product remainder: ", bound_summary(remainder(a * b)))
    if profiling
        println("\nDense product CPU profile:")
        profile_workload(() -> a * b)
    end

    x = taylor_models(
        fill(interval(Float64, -1 // 50, 1 // 50), 6);
        order = 3, table_bytes
    )
    println("Uncertainty expression: ", measure_kernel(() -> uncertainty_map(x)))
    model = uncertainty_map(x)
    println(
        "Uncertainty range / remainder: ",
        (bound_summary(enclose(model)), bound_summary(remainder(model)))
    )
    if profiling
        println("\nUncertainty expression CPU profile:")
        profile_workload(() -> uncertainty_map(x))
        Profile.Allocs.clear()
        Profile.Allocs.@profile sample_rate = 1.0 uncertainty_map(x)
        allocations = Profile.Allocs.fetch().allocs
        owners = Dict{String, Tuple{Int, Int}}()
        for allocation in allocations
            index = findfirst(
                frame -> occursin("interval_models.jl", string(frame.file)) ||
                    occursin("arithmetic.jl", string(frame.file)), allocation.stacktrace
            )
            index === nothing && continue
            owner = string(allocation.stacktrace[index])
            count, bytes = get(owners, owner, (0, 0))
            owners[owner] = (count + 1, bytes + allocation.size)
        end
        println("\nUncertainty allocations by first model/arithmetic frame:")
        for pair in sort!(collect(owners); by = pair -> last(pair)[2], rev = true)
            println(pair)
        end
    end
    return nothing
end

for budget in (0, 32 * 1024^2)
    report(; table_bytes = budget, profiling = "--profile" in ARGS)
end
