# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Saving and restoring an evaluation map
#
# Compile polynomials into a reusable numeric map, then use Julia's standard
# Serialization library for a local checkpoint. This is Julia-specific storage,
# not a portable interchange format or a long-term archival guarantee.
# Only deserialize files from a trusted source.
using DifferentialAlgebra
using Serialization
using CairoMakie

x1, x2, x3, x4 = variables((:x1, :x2, :x3, :x4); order = 5)
p = 2 - 0.75x1 + 1.2x2^2 + 0.1x1 * x4 - 0.05x3^3
polynomials = [1 + x1 + 0.5x2^2 - 0.2x3 * x4, -2 + 0.3x2 + x1 * x3 + 0.8x4^2]
snapshot = (; scalar = compile(p), map = compile(polynomials))
points = [[0.0, 0.0, 0.0, 0.0], [0.1, -0.2, 0.3, -0.4], [-1.0, 0.5, 0.25, 1.0]]
expected = [(scalar = [p(x)], map = evaluate(polynomials, x)) for x in points]
println("Scalar polynomial: ", p)
println("Vector map: ", polynomials)

# A temporary directory keeps running this example free of persistent artifacts.
# Numeric evaluation of a compiled map survives a new algebra initialization.
loaded = mktempdir() do directory
    path = joinpath(directory, "map.jls")
    serialize(path, snapshot)
    deserialize(path)
end
variables(1; order = 1)
for (point, reference) in zip(points, expected)
    @assert loaded.scalar(point) == reference.scalar
    @assert loaded.map(point) == reference.map
    println((point = point, scalar = only(loaded.scalar(point)), map = loaded.map(point)))
end
println("Scalar and vector checkpoints agree at all validation points.")

# ## The restored map preserves the entire stored polynomial
# Follow a curve through the four-dimensional input space. Both compiled
# snapshots remain evaluable after the global algebra has been replaced.
# This checks checkpoint equality, not error relative to a generating function.
parameter = range(-1, 1; length = 201)
path = [[t, 0.3sin(2t), 0.2cos(t), -0.4t] for t in parameter]
before_scalar = [only(snapshot.scalar(v)) for v in path]
after_scalar = [only(loaded.scalar(v)) for v in path]
before_map, after_map = snapshot.map.(path), loaded.map.(path)
@assert before_scalar == after_scalar && before_map == after_map
fig = Figure(size = (1050, 420), fontsize = 14)
ax = Axis(fig[1, 1]; xlabel = "input-curve parameter", ylabel = "scalar polynomial", title = "Checkpoint and algebra reinitialization")
lines!(ax, parameter, before_scalar; color = :black, linewidth = 2, label = "Original compiled snapshot")
scatter!(ax, parameter[1:10:end], after_scalar[1:10:end]; color = :dodgerblue, markersize = 7, label = "Restored snapshot")
axislegend(ax; position = :rt, labelsize = 11)
ax = Axis(fig[1, 2]; xlabel = "map component 1", ylabel = "map component 2", title = "Image of the same 4D input curve")
lines!(ax, first.(before_map), last.(before_map); color = :black, linewidth = 2)
scatter!(ax, first.(after_map[1:10:end]), last.(after_map[1:10:end]); color = :dodgerblue, markersize = 7)
fig
