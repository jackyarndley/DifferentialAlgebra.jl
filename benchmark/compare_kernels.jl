# Compare arithmetic and analytic Kepler propagation with an older checkout in
# the same process. Instantiate both checkouts first. No ODE extensions are loaded
# for the older module, so this deliberately does not compare SciML integration.
# Run `julia --project=benchmark benchmark/compare_kernels.jl OLD_CHECKOUT output.toml`.
include("orbits.jl")
using .OrbitBenchmarks, TOML

length(ARGS) == 2 || error("Provide an older checkout and an output TOML path")
const BASELINE = abspath(first(ARGS))
push!(LOAD_PATH, BASELINE)
module Before
    include(joinpath(Main.BASELINE, "src", "DifferentialAlgebra.jl"))
end
const B = OrbitBenchmarks
const Old = Before.DifferentialAlgebra
B.primal(p::Old.TaylorPolynomial) = Old.constant_term(p)
B.polynomial_coefficient(p::Old.TaylorPolynomial, powers) = Old.coefficient(p, powers)

function compare()
    data = B.metadata(11)
    data["baseline_src_sha256"] = bytes2hex(B.sha256(join(read(joinpath(BASELINE, "src", f), String) for f in sort(readdir(joinpath(BASELINE, "src"))))))
    function record(name, order, functions; scale)
        old, new = functions[1](), functions[2]()
        if old isa Real
            old, new = (old,), (new,)
        end
        errors = B.coefficient_errors(old, new, order; perturbation_scale = scale)
        @assert maximum(errors) < 1.0e-10 errors
        timings = B.measure(functions; samples = 11)
        println(order, " ", name, ": ", timings[1]["seconds"] / timings[2]["seconds"], "×")
        flush(stdout)
        return push!(data["cases"], Dict("order" => order, "method" => name, "before" => timings[1], "after" => timings[2], "coefficient_errors" => errors))
    end
    for order in (4, 8, 12)
        a0 = exp(sum(Old.variables(6; order)) / 10); b0 = 2 + a0 * a0
        a1 = exp(sum(B.DA.variables(6; order)) / 10); b1 = 2 + a1 * a1
        for (name, f) in (
                ("quotient", (a, b) -> a / b), ("multiply", (a, b) -> a * b),
                ("sqrt", (a, b) -> sqrt(a)), ("weighted_sum", (a, b) -> muladd(0.1, a, b)),
            )
            record(name, order, (() -> f(a0, b0), () -> f(a1, b1)); scale = 0.1)
        end
        for s in B.scenarios()[2:3]
            u0 = s.u0 + B.PERTURBATION_SCALE * B.SVector{6}(Old.variables(6; order))
            u1 = B.initial(s.u0, order, :DifferentialAlgebra)
            functions = (() -> B.kepler(u0, s.tf, order), () -> B.kepler(u1, s.tf, order))
            record(s.name, order, functions; scale = B.PERTURBATION_SCALE)
        end
    end
    return open(last(ARGS), "w") do io
        TOML.print(io, data)
    end
end
compare()
