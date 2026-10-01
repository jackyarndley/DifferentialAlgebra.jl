# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex14_ads_gaussian.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Automatic domain splitting of a Gaussian
#
# ADS builds local Taylor maps on a physical box. It refines patches whose
# estimated truncation error exceeds the tolerance. A separate grid tests
# actual error; a successful heuristic estimate is not a certified bound.
using DifferentialAlgebra
using CairoMakie

gaussian(x) = exp(-(x[1]^2 + x[2]^2) / 2)
tolerance = 1.0e-4
map = adaptive_map(gaussian, [-3.0, -3.0], [3.0, 3.0]; order = 5, atol = tolerance)
grid = range(-3, 3; length = 101)
exact = [gaussian([x, y]) for x in grid, y in grid]
approximation = [map([x, y]) for x in grid, y in grid]
errors = abs.(approximation - exact)
@assert map.converged && length(map.patches) > 1
@assert maximum(errors) < tolerance
println(map)
println((maximum_error = maximum(errors), rms_error = sqrt(sum(abs2, errors) / length(errors))))

# The partition concentrates small patches where the map needs higher resolution.
fig = Figure(size = (1150, 440), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", title = "Gaussian map", aspect = DataAspect())
values = heatmap!(ax, grid, grid, approximation; colormap = :viridis)
Colorbar(fig[2, 1], values; vertical = false, label = "Gaussian value")
partition = Axis(fig[1, 2]; xlabel = "x", ylabel = "y", title = "$(length(map.patches)) Taylor patches", aspect = DataAspect())
for patch in map.patches
    x0, y0 = patch.lower
    x1, y1 = patch.upper
    lines!(partition, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color = :steelblue, linewidth = 1)
end
err = Axis(fig[1, 3]; xlabel = "x", ylabel = "y", title = "Independent grid error", aspect = DataAspect())
heat = heatmap!(err, grid, grid, log10.(max.(errors, 1.0e-12)); colormap = :magma, colorrange = (-12, -4))
Colorbar(fig[1, 4], heat; label = "log₁₀ absolute error")
fig
