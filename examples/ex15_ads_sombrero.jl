# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex15_ads_sombrero.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # ADS across the sombrero's origin
#
# The smooth function sin(√q)/√q has a removable singularity at q=0.
# Adding a small number under the square root changes the function and
# introduces bias. Instead, use its entire series for a patch centered at zero.
using DifferentialAlgebra
using CairoMakie

function sombrero(x)
    q = x[1]^2 + x[2]^2
    if iszero(constant_term(q))
        q isa TaylorPolynomial || return one(q)
        value = term = one(q)
        for k in 1:truncation_order()
            term *= -q / ((2k) * (2k + 1))
            value += term
        end
        return value
    end
    r = sqrt(q)
    return sin(r) / r
end

tolerance = 1.0e-5
map = adaptive_map(sombrero, [-6.0, -6.0], [6.0, 6.0]; order = 7, atol = tolerance)
grid = range(-6, 6; length = 121)
reference = [sinc(hypot(x, y) / π) for x in grid, y in grid]
values = [map([x, y]) for x in grid, y in grid]
errors = abs.(values - reference)
@assert map.converged && length(map.patches) > 1
@assert abs(map([0.0, 0.0]) - 1) < tolerance
@assert maximum(errors) < tolerance
println(map)
println((origin = map([0.0, 0.0]), maximum_error = maximum(errors)))

#-
fig = Figure(size = (1100, 470), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", title = "Sombrero and ADS partition", aspect = DataAspect())
surface = heatmap!(ax, grid, grid, values; colormap = :viridis)
Colorbar(fig[2, 1], surface; vertical = false, label = "sin(r) / r")
for patch in map.patches
    x0, y0 = patch.lower
    x1, y1 = patch.upper
    lines!(ax, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color = (:white, 0.6), linewidth = 0.7)
end
err = Axis(fig[1, 2]; xlabel = "x", ylabel = "y", title = "Error against unregularized sinc", aspect = DataAspect())
heat = heatmap!(err, grid, grid, log10.(max.(errors, 1.0e-13)); colormap = :magma, colorrange = (-13, -5))
Colorbar(fig[1, 3], heat; label = "log₁₀ absolute error")
fig
