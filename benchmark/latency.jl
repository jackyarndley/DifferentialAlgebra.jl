# Run in fresh processes, with the benchmark environment instantiated:
# julia --startup-file=no --project=benchmark benchmark/latency.jl results.toml
# JULIA_DEPOT_PATH and PrecompileTools preferences are inherited by children.
using Dates, SHA, TOML

function latency_workload(point)
    x, y = variables(2; order = 5)
    p = exp(x / 10) * sin(y / 10) + sqrt(2 + x - y) / (3 + x + y)
    map = CompiledMap([p, differentiate(p, 1)])
    return evaluate(map, point)
end

function latency_build()
    built = @timed Base.compilecache(Base.identify_package("DifferentialAlgebra"))
    paths = filter(path -> path isa String && isfile(path), collect(built.value))
    return Dict(
        "build_seconds" => built.time, "cache_bytes" => sum(filesize, paths),
        "cache_files" => paths
    )
end

function latency_evaluate(load_seconds)
    @assert !DifferentialAlgebra.isinitialized()
    point = [0.1, 0.2]
    # Keep workload compilation inside the timer, rather than compiling it as
    # part of this measurement function before the timer starts.
    first_call = @timed Base.invokelatest(latency_workload, point)
    warm_call = @timed Base.invokelatest(latency_workload, point)
    @assert first_call.value ≈ warm_call.value
    return Dict(
        "load_seconds" => load_seconds,
        "first_seconds" => first_call.time, "first_compile_seconds" => first_call.compile_time,
        "warm_seconds" => warm_call.time, "first_bytes" => first_call.bytes
    )
end

function latency_main(output; samples = 3)
    root = dirname(@__DIR__)
    files = sort(filter(f -> endswith(f, ".jl"), readdir(joinpath(root, "src"); join = true)))
    source_hash = bytes2hex(sha256(join(read.(files, String))))
    results = Dict{String, Any}(
        "julia" => string(VERSION), "date" => string(now(UTC)),
        "os" => string(Sys.KERNEL), "cpu" => Sys.cpu_info()[1].model,
        "source_sha256" => source_hash, "samples" => samples
    )
    project = dirname(Base.active_project())
    preference_file = joinpath(project, "LocalPreferences.toml")
    preferences = isfile(preference_file) ? TOML.parsefile(preference_file) : Dict()
    results["precompile_workload"] = get(get(preferences, "DifferentialAlgebra", Dict()), "precompile_workload", true)
    for mode in ("build", "core")
        records = Dict[]
        for _ in 1:samples
            mktemp() do path, io
                close(io)
                run(`$(Base.julia_cmd()) --startup-file=no --project=$project $(@__FILE__) --child $mode $path`)
                push!(records, TOML.parsefile(path))
            end
        end
        results[mode] = records
    end
    open(io -> TOML.print(io, results), output, "w")
    return println("Latency results: ", abspath(output))
end

if !isempty(ARGS) && first(ARGS) == "--child"
    if ARGS[2] == "build"
        result = latency_build()
    else
        loaded = @timed using DifferentialAlgebra
        result = latency_evaluate(loaded.time)
        result["load_compile_seconds"] = loaded.compile_time
    end
    open(io -> TOML.print(io, result), ARGS[3], "w")
else
    latency_main(isempty(ARGS) ? "latency.toml" : first(ARGS))
end
