# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Kepler's equation by coupled map inversion
#
# Solve M = E - e sin(E) near M=π/2 and e=0.3. Inverting the augmented
# map (δE, δe) ↦ (E - e sin(E) - M, δe) gives the solution's dependence
# on both mean anomaly and eccentricity. Check its sensitivities analytically,
# then compare a fixed-M slice with independently converged scalar roots.
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
joint_solution = E0 + inverse[1]
solution = E0 + evaluate(inverse[1], [zero(δe), δe])
residual = coefficient_norm(solution - e * sin(solution) - mean)
@assert residual < 1.0e-10
println("E(δe) = ", solution)
println("Maximum implicit-equation coefficient residual: ", residual)

# ## Sensitivities without solving again
#
# The inverse's first coordinate now means δM, not the original δE.
# Implicit differentiation gives ∂E/∂M = 1/(1-e cos(E)) and
# ∂E/∂e = sin(E)/(1-e cos(E)) at the nominal root.
sensitivity = constant_term(gradient(joint_solution))
exact_sensitivity = [1.0, sin(E0)] / (1 - e0 * cos(E0))
@assert maximum(abs, sensitivity - exact_sensitivity) < 1.0e-13
println("Sensitivity to (M, e), DA / analytical: ", (sensitivity, exact_sensitivity))

# ## A finite uncertainty interval
#
# Setting δM=0 yields E(δe). Scalar root solves provide an independent check
# throughout e ∈ [0.1,0.5]; this sampled discrepancy is not a uniform bound.
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
