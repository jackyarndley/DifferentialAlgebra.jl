# Sample the warmed polynomial propagation paths on the current implementation.
# Run `julia --project=benchmark benchmark/profile_orbits.jl [output.txt]`.
# Compilation and global-basis construction are excluded from the CPU profile.
using Profile
include("orbits.jl")
using .OrbitBenchmarks

function profile_orbits(io)
    B = OrbitBenchmarks
    scenario = B.scenarios()[2]
    for (order, analytic) in ((1, false), (8, false), (12, true))
        u0 = B.initial(scenario.u0, order, :DifferentialAlgebra)
        f = analytic ? (() -> B.kepler(u0, scenario.tf, order)) :
            (() -> B.propagate(u0, scenario.tf, scenario.j2, scenario.algorithm, scenario.tolerance))
        f(); f()
        GC.gc()
        Profile.clear()
        Profile.@profile for _ in 1:(order == 1 ? 100 : 5)
            f()
        end
        println(io, "\n", analytic ? "Kepler" : "Vern9", ", order ", order)
        Profile.print(IOContext(io, :displaysize => (200, 220)); format = :flat, sortedby = :count, mincount = 5, threads = Threads.threadid())
    end
    return nothing
end

if isempty(ARGS)
    profile_orbits(stdout)
else
    open(profile_orbits, only(ARGS), "w")
end
