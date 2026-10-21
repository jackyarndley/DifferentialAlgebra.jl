# julia --startup-file=no --project=benchmark benchmark/interval_kernels.jl --profile
using DifferentialAlgebra, IntervalArithmetic
using Profile

include("measurement.jl")

const DA = DifferentialAlgebra
const IA = IntervalArithmetic
const EXT = Base.get_extension(DA, :DifferentialAlgebraIntervalArithmeticExt)

uncertainty_map(x) = exp(x[1] + x[2] / 2) * cos(x[3] - x[4]) +
    log(2 + x[5] * x[6])

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
    powers = EXT.power_cache(a._coordinates.normalized, 6, Interval{Float64})
    println("Dense model product: ", measure(() -> a * b))
    println(
        "Dense polynomial range: ", measure(
            () -> EXT.polynomial_bound(
                p, powers,
                Interval{Float64}
            )
        )
    )
    println(
        "Discarded products: ", measure(
            () -> EXT.discarded_bound(
                p, q, 3,
                powers, Interval{Float64}
            )
        )
    )
    println("Dense model range: ", enclose(a * b))
    if profiling
        println("\nDense product CPU profile:")
        profile_workload(() -> a * b)
    end

    x = taylor_models(
        fill(interval(Float64, -1 // 50, 1 // 50), 6);
        order = 3, table_bytes
    )
    println("Uncertainty expression: ", measure(() -> uncertainty_map(x)))
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
