"""
Adaptive orbit and high-order Kepler benchmarks with independent accuracy checks.
Run `julia --project=benchmark benchmark/orbits.jl [--quick] output.toml`.
All backends use the same model, initial state, algorithm and nominal error norm.
"""
module OrbitBenchmarks
import DifferentialAlgebra as DA
import TaylorSeries as TS
import DifferentiationInterface as DI
import ForwardDiff
using LinearAlgebra, StaticArrays, Statistics, Printf, TOML, SHA, Dates
using BenchmarkTools
using OrdinaryDiffEqTsit5: OrdinaryDiffEqTsit5, Tsit5
using OrdinaryDiffEqVerner: OrdinaryDiffEqVerner, Vern9
using SciMLBase: SciMLBase, ODEProblem, solve

include("orbit_models.jl")
const PERTURBATION_SCALE = 1.0e-3

function initial(u0, order, backend)
    x = backend === :DifferentialAlgebra ? DA.variables(6; order) :
        TS.variables!(Float64, [Symbol("x", i) for i in 1:6]; order, nowarn = true)
    return u0 + PERTURBATION_SCALE * SVector{6}(x)
end

function state_stm(polynomials::AbstractVector{<:DA.TaylorPolynomial})
    return DA.constant_term(polynomials), DA.linear_part(polynomials) / PERTURBATION_SCALE
end
function state_stm(polynomials::AbstractVector{<:TS.TaylorN})
    J = zeros(6, 6)
    for j in 1:6
        powers = zeros(Int, 6); powers[j] = 1
        for i in 1:6
            J[i, j] = TS.getcoeff(polynomials[i], powers) / PERTURBATION_SCALE
        end
    end
    return TS.constant_term.(polynomials), J
end

polynomial_coefficient(p::DA.TaylorPolynomial, powers) = DA.coefficient(p, powers)
polynomial_coefficient(p::TS.TaylorN, powers) = TS.getcoeff(p, Int.(powers))
function coefficient_errors(a, b, order; perturbation_scale = PERTURBATION_SCALE)
    @assert length(a) == length(b)
    errors, magnitudes = zeros(order + 1), zeros(order + 1)
    for powers in DA.multiindices(order, 6)
        d = Int(sum(powers)); scale = perturbation_scale^d
        for i in eachindex(a)
            av = polynomial_coefficient(a[i], powers) / scale
            bv = polynomial_coefficient(b[i], powers) / scale
            errors[d + 1] = max(errors[d + 1], abs(av - bv))
            magnitudes[d + 1] = max(magnitudes[d + 1], abs(av), abs(bv))
        end
    end
    return errors ./ max.(1.0, magnitudes)
end

relative_error(a, b) = maximum(abs, a - b) / max(1.0, maximum(abs, b))
function validate_stm(value, reference; tolerance)
    state_error, stm_error = relative_error(value[1], reference[1]), relative_error(value[2], reference[2])
    @assert state_error < tolerance (state_error, tolerance)
    @assert stm_error < tolerance (stm_error, tolerance)
    return (state_error, stm_error)
end

# Compilation, basis construction and AD preparation are excluded;
# problem/solver creation and extraction of the final STM are included.
# Rotate backend order, use warmed batched samples for fast cases, and retain GC.
function measure(functions; samples)
    times = [Float64[] for _ in functions]
    bytes = [Int[] for _ in functions]
    allocations = [Int[] for _ in functions]
    batches = Int[]
    benchmarks = map(f -> @benchmarkable($f()), functions)
    for f in functions
        f(); f()
        GC.gc()
        trial = @timed f()
        push!(batches, clamp(ceil(Int, 0.005 / max(trial.time, 1.0e-9)), 1, 1000))
    end
    for sample in 1:samples, shift in 0:(length(functions) - 1)
        j = mod1(sample + shift, length(functions))
        GC.gc()
        trial = BenchmarkTools.run(benchmarks[j]; samples = 1, evals = batches[j], seconds = 3600, gctrial = false, gcsample = false)
        estimate = median(trial)
        push!(times[j], estimate.time / 1.0e9)
        push!(bytes[j], estimate.memory)
        push!(allocations[j], estimate.allocs)
    end
    return [
        Dict(
            "seconds" => median(times[i]), "bytes" => round(Int, median(bytes[i])),
            "allocations" => round(Int, median(allocations[i])), "samples" => times[i], "batch" => batches[i]
        ) for i in eachindex(functions)
    ]
end

function row(scenario, method, backend, order, timing, errors; solution = nothing, coefficient_error = Float64[])
    result = merge(
        timing, Dict(
            "scenario" => scenario.name, "method" => method,
            "backend" => backend, "order" => order, "duration" => scenario.tf,
            "tolerance" => scenario.tolerance, "state_error" => errors[1],
            "coefficient_errors" => coefficient_error
        )
    )
    errors[2] === nothing || (result["stm_error"] = errors[2])
    if solution !== nothing
        merge!(
            result, Dict(
                "accepted_steps" => solution.stats.naccept,
                "rejected_steps" => solution.stats.nreject, "rhs_calls" => solution.stats.nf
            )
        )
    end
    @printf(
        "%-20s %-7s %-20s order=%2d %10.3f ms %12d B state=%.2g STM=%s\n",
        scenario.name, method, backend, order, 1000timing["seconds"], timing["bytes"], errors[1],
        errors[2] === nothing ? "n/a" : @sprintf("%.2g", errors[2])
    )
    flush(stdout)
    return result
end

function run_scenario(scenario; orders, samples, analytic = false)
    s = scenario
    method = analytic ? "Kepler" : string(nameof(typeof(s.algorithm)))
    flow = analytic ? (u -> kepler(u, s.tf, 1)) : (u -> propagate(u, s.tf, s.j2, s.algorithm, s.tolerance).u[end])
    reference = variational_reference(s.u0, s.tf, s.j2)
    if iszero(s.j2)
        # An independently solved Kepler equation also verifies the reference.
        exact = DI.value_and_jacobian(u -> kepler(u, s.tf, 1), DI.AutoForwardDiff(chunksize = 6), s.u0)
        validate_stm(reference, exact; tolerance = 2.0e-9)
        reference = exact
    end
    backend = DI.AutoForwardDiff(chunksize = 6)
    prep = DI.prepare_jacobian(flow, backend, s.u0)
    numeric = () -> (flow(s.u0), nothing)
    differentiated = () -> DI.value_and_jacobian(flow, prep, backend, s.u0)
    error_limit = 1000s.tolerance
    ad_errors = validate_stm(differentiated(), reference; tolerance = analytic ? 1.0e-12 : error_limit)
    numeric_error = relative_error(numeric()[1], reference[1])
    @assert numeric_error < (analytic ? 1.0e-12 : error_limit)
    timings = measure((numeric, differentiated); samples)
    scalar_solution = analytic ? nothing : propagate(s.u0, s.tf, s.j2, s.algorithm, s.tolerance)
    # Collect solver statistics with an equivalent fully seeded dual state;
    # the timed AD call above always goes through DifferentiationInterface.
    dual_input = SVector{6}(ntuple(i -> ForwardDiff.Dual(s.u0[i], ntuple(j -> Float64(i == j), 6)), 6))
    dual_solution = analytic ? nothing : propagate(dual_input, s.tf, s.j2, s.algorithm, s.tolerance)
    rows = [
        row(s, method, "Float64", 0, timings[1], (numeric_error, nothing); solution = scalar_solution),
        row(s, method, "DI-ForwardDiff", 1, timings[2], ad_errors; solution = dual_solution),
    ]
    for order in orders
        native, taylor = initial(s.u0, order, :DifferentialAlgebra), initial(s.u0, order, :TaylorSeries)
        propagate_native = analytic ? (() -> kepler(native, s.tf, order)) : (() -> propagate(native, s.tf, s.j2, s.algorithm, s.tolerance).u[end])
        propagate_taylor = analytic ? (() -> kepler(taylor, s.tf, order)) : (() -> propagate(taylor, s.tf, s.j2, s.algorithm, s.tolerance).u[end])
        a, b = propagate_native(), propagate_taylor()
        errors = coefficient_errors(a, b, order)
        # Roundoff can change the adaptive mesh even when step counts agree.
        @assert maximum(errors) < (analytic ? 1.0e-10 : 1.0e-8) (scenario = s.name, order, errors)
        ae, be = validate_stm(state_stm(a), reference; tolerance = error_limit), validate_stm(state_stm(b), reference; tolerance = error_limit)
        # Test a nonzero perturbation against scalar propagation, not just nominal values.
        delta = SVector(0.02, -0.015, 0.01, -0.02, 0.01, 0.015)
        perturbed = flow(s.u0 + PERTURBATION_SCALE * delta)
        perturbed_errors = (relative_error(DA.evaluate(a, delta), perturbed), relative_error([TS.evaluate(p, delta) for p in b], perturbed))
        @assert maximum(perturbed_errors) < (order == 1 ? 5.0e-4 : 5.0e-7) perturbed_errors
        if !analytic
            reference_map = iszero(s.j2) ? kepler(native, s.tf, order) : propagate(native, s.tf, s.j2, Vern9(), s.tolerance / 20).u[end]
            flow_errors = (coefficient_errors(a, reference_map, order), coefficient_errors(reference_map, b, order))
            @assert max(maximum(flow_errors[1]), maximum(flow_errors[2])) < 5.0e-7 (scenario = s.name, order, flow_errors)
        else
            flow_errors = (Float64[], Float64[])
        end
        # First order includes final-state/STM extraction for an equal AD comparison.
        fd = order == 1 ? (() -> state_stm(propagate_native())) : (() -> (propagate_native(), nothing))
        ft = order == 1 ? (() -> state_stm(propagate_taylor())) : (() -> (propagate_taylor(), nothing))
        timings = measure((fd, ft); samples)
        da_solution = analytic ? nothing : propagate(native, s.tf, s.j2, s.algorithm, s.tolerance)
        ts_solution = analytic ? nothing : propagate(taylor, s.tf, s.j2, s.algorithm, s.tolerance)
        for (i, (name, timing, accuracy, solution)) in enumerate((("DifferentialAlgebra", timings[1], ae, da_solution), ("TaylorSeries", timings[2], be, ts_solution)))
            result = row(s, method, name, order, timing, accuracy; solution, coefficient_error = errors)
            result["reference_coefficient_errors"] = flow_errors[i]
            result["perturbed_error"] = perturbed_errors[i]
            push!(rows, result)
        end
    end
    return rows
end

function metadata(samples)
    root = dirname(@__DIR__)
    sources = sort(vcat(readdir(joinpath(root, "src"); join = true), readdir(joinpath(root, "ext"); join = true)))
    source_hash = bytes2hex(sha256(join(read(file, String) for file in sources)))
    return Dict(
        "core_source_sha256" => source_hash,
        "utc" => string(now(UTC)), "cpu_model" => first(Sys.cpu_info()).model,
        "julia" => string(VERSION), "cpu" => Sys.CPU_NAME, "kernel" => string(Sys.KERNEL),
        "threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "samples" => samples, "perturbation_scale" => PERTURBATION_SCALE,
        "benchmarktools" => string(pkgversion(BenchmarkTools)),
        "differentialalgebra" => string(pkgversion(DA)), "taylorseries" => string(pkgversion(TS)),
        "differentiationinterface" => string(pkgversion(DI)), "forwarddiff" => string(pkgversion(ForwardDiff)),
        "ordinarydiffeqtsit5" => string(pkgversion(OrdinaryDiffEqTsit5)),
        "ordinarydiffeqverner" => string(pkgversion(OrdinaryDiffEqVerner)),
        "scimlbase" => string(pkgversion(SciMLBase)), "cases" => Dict[]
    )
end

function main()
    quick = "--quick" in ARGS
    samples = quick ? 3 : 7
    data = metadata(samples)
    println("Orbit benchmarks: Julia ", VERSION, ", ", Sys.CPU_NAME, ", ", Sys.KERNEL)
    output = isempty(ARGS) || !endswith(last(ARGS), ".toml") ? nothing : last(ARGS)
    for analytic in (false, true), (i, scenario) in enumerate(scenarios())
        analytic && i ∉ (2, 3) && continue
        orders = analytic ? (quick ? (1, 4, 8, 12) : (1, 2, 4, 6, 8, 10, 12)) : (quick ? (1, 4, 8) : (1, 2, 4, 6, 8))
        append!(data["cases"], run_scenario(scenario; orders, samples, analytic))
        if output !== nothing
            open(output, "w") do io
                TOML.print(io, data)
            end
        end
    end
    return data
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
