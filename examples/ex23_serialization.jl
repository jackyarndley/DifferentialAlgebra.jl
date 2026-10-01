# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex23_dace_io_roundtrip.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Saving and restoring an evaluation map
#
# Compile polynomials into a reusable numeric map, then use Julia's standard
# Serialization library for a local checkpoint. This is Julia-specific storage,
# not a portable interchange format or a long-term archival guarantee.
# Only deserialize files from a trusted source.
using DifferentialAlgebra
using Serialization

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
