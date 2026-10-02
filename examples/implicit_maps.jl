# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Kepler's equation by coupled map inversion
#
# Solve M = E - e sin(E), with M=π/2 and uncertain eccentricity near 0.3.
# First locate the scalar root, then invert the augmented map
# (δE, δe) ↦ (E - e sin(E) - M, δe). Setting its first output to zero
# gives E as a polynomial in δe. The extra coordinate makes the map square.
using DifferentialAlgebra
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

function kepler_root(mean, eccentricity)
    E = mean
    for _ in 1:40
        step = (E - eccentricity * sin(E) - mean) / (1 - eccentricity * cos(E))
        E -= step
        abs(step) < 2.0e-15 && return E
    end
    error("Kepler root did not converge")
end

δE, δe = variables((:δE, :δe); order = 20)
mean, e0 = π / 2, 0.3
E0 = kepler_root(mean, e0)
E, e = E0 + δE, e0 + δe
constraint = E - e * sin(E) - mean
inverse = invert([constraint, δe])
solution = E0 + evaluate(inverse[1], [zero(δe), δe])
residual = coefficient_norm(solution - e * sin(solution) - mean)
@assert residual < 1.0e-10
println("E(δe) = ", solution)
println("Maximum implicit-equation coefficient residual: ", residual)

# Compare the polynomial with independently converged scalar roots.
eccentricities = range(0.1, 0.5; length = 151)
reference = [kepler_root(mean, e) for e in eccentricities]
values = [solution([0.0, e - e0]) for e in eccentricities]
@assert maximum(abs, values - reference) < 1.0e-10
println("Maximum error on e ∈ [0.1, 0.5]: ", maximum(abs, values - reference))

#-
fig = Figure(size = (1000, 390), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "Eccentricity", ylabel = "Eccentric anomaly (rad)", title = "Implicit Kepler map")
lines!(ax, eccentricities, reference; color = :black, linewidth = 3, label = "Scalar Newton")
lines!(ax, eccentricities, values; linestyle = :dash, linewidth = 2, label = "Order 20 inverse map")
axislegend(ax; position = :lt)
err = Axis(fig[1, 2]; xlabel = "Eccentricity", ylabel = "Absolute error (rad)", yscale = log10)
lines!(err, eccentricities, max.(abs.(values - reference), eps(Float64)); linewidth = 2)
fig
