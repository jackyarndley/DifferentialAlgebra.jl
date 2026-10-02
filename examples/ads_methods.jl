# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Choosing an automatic domain splitting method
#
# A Gaussian on an anisotropic box illustrates how the error estimator and
# split direction affect the partition. All methods use the same interface,
# polynomial order, tolerance and independent validation grid.
using DifferentialAlgebra
using CairoMakie

gaussian(x) = exp(-(4x[1]^2 + x[2]^2) / 2)
lower, upper = [-2.0, -2.0], [2.0, 2.0]
tolerance = 1.0e-4
methods = [
    ("Guard degrees", GuardedTail(), :tail),
    ("Coefficient decay", ExtrapolatedTail(), :tail),
    ("Last two degrees", LastTerms(), :tail),
    ("Longest relative side", GuardedTail(), :width),
]
maps = [
    adaptive_map(gaussian, lower, upper; order = 5, atol = tolerance, estimator, splitter)
        for (_, estimator, splitter) in methods
]

# GuardedTail computes two extra degrees. ExtrapolatedTail fits an exponential
# decay to retained coefficient maxima, and LastTerms inspects the last two
# degrees without fitting. The :tail splitter favors directions that reduce
# large coefficients; :width provides a simple geometric alternative.
# These are heuristic indicators. Measure error independently of their probes.
grid = range(-2, 2; length = 101)
reference = [gaussian([x, y]) for x in grid, y in grid]
errors = [abs.([map([x, y]) for x in grid, y in grid] - reference) for map in maps]
for ((name, _, _), map, error) in zip(methods, maps, errors)
    @assert map.converged
    @assert maximum(error) < tolerance
    println((method = name, patches = length(map.patches), maximum_error = maximum(error)))
end

# Existing partitions can be refined without merging their boundaries.
# Always provide the original function so new expansions recover missing terms.
refined = adaptive_map(gaussian, maps[1]; atol = tolerance / 10)
refined_error = maximum(abs(refined([x, y]) - gaussian([x, y])) for x in grid, y in grid)
@assert refined.converged && refined_error < tolerance / 10
println((refined_patches = length(refined.patches), maximum_error = refined_error))

# The upper row shows partitions in physical coordinates. The lower row uses
# the same color scale for independently measured errors in all four methods.
fig = Figure(size = (1300, 710), fontsize = 14)
for (column, ((name, _, _), map, error)) in enumerate(zip(methods, maps, errors))
    domain = Axis(
        fig[1, column]; xlabel = "x", ylabel = "y",
        title = "$name\n$(length(map.patches)) patches", aspect = DataAspect()
    )
    for patch in map.patches
        x0, y0 = patch.lower
        x1, y1 = patch.upper
        lines!(domain, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color = :steelblue, linewidth = 0.7)
    end
    accuracy = Axis(fig[2, column]; xlabel = "x", ylabel = "y", title = "Independent grid error", aspect = DataAspect())
    heatmap!(accuracy, grid, grid, log10.(max.(error, 1.0e-12)); colormap = :magma, colorrange = (-12, -4))
end
Colorbar(fig[3, 1:4]; colormap = :magma, limits = (-12, -4), vertical = false, label = "log₁₀ absolute error")
fig
