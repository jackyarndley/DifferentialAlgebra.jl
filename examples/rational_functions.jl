# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # A rational function near a nonzero point
#
# The variable δx represents a displacement, not the physical coordinate.
# For f(x) = x/(x²+1), f(3) = 0.3 and f′(3) = -0.08.
using DifferentialAlgebra
using CairoMakie

δx, = variables((:δx,); order = 1)
x = 3 + δx
f = 1 / (x + 1 / x)
@assert constant_term(f) ≈ 0.3
@assert coefficient(f, [1]) ≈ -0.08
println("x = ", x)
println("f(3 + δx) = ", f)
println("Linear prediction at x = 3.01: ", f(0.01))
println("Direct value: ", 3.01 / (3.01^2 + 1))

# ## Displacement coordinates and the useful local range
# The horizontal coordinate below is physical x; the polynomial argument is
# always x-3. Higher order improves the local fit, but does not remove the
# complex poles at ±i. The Taylor series about 3 has convergence radius √10.
δx, = variables((:δx,); order = 12)
orders = (1, 3, 6, 12)
polynomials = [with_order(() -> (3 + δx) / (1 + (3 + δx)^2), n) for n in orders]
physical = range(1, 5; length = 301)
reference = physical ./ (1 .+ physical .^ 2)
values = [[p(x - 3) for x in physical] for p in polynomials]
fig = Figure(size = (1050, 400), fontsize = 14)
ax = Axis(fig[1, 1]; xlabel = "physical x", ylabel = "x/(1+x²)", title = "Local expansion about x=3")
lines!(ax, physical, reference; color = :black, linewidth = 3, label = "Original function")
err = Axis(fig[1, 2]; xlabel = "physical x", ylabel = "sampled absolute error", yscale = log10, title = "Order and distance from the center")
for (n, v) in zip(orders, values)
    lines!(ax, physical, v; label = "Order $n")
    lines!(err, physical, max.(abs.(v - reference), eps(Float64)); label = "Order $n")
end
vlines!(ax, [3]; color = :gray, linestyle = :dot)
vlines!(err, [3]; color = :gray, linestyle = :dot)
axislegend(ax; position = :rt, labelsize = 11)
fig

# The next plot crosses the convergence radius. A finite polynomial remains
# evaluable there, but increasing its order need not improve the approximation.
# Error floors in log plots are display conventions, not uniform bounds.
displacements = range(0, 5; length = 301)
radius = sqrt(10.0)
convergence_fig = Figure(size = (720, 400), fontsize = 14)
ax = Axis(convergence_fig[1, 1]; xlabel = "displacement δx (physical x=3+δx)", ylabel = "sampled absolute error", yscale = log10, title = "Beyond the Taylor convergence radius")
for (n, p) in zip(orders, polynomials)
    errors = [abs(p(t) - (3 + t) / (1 + (3 + t)^2)) for t in displacements]
    lines!(ax, displacements, max.(errors, eps(Float64)); label = "Order $n")
end
vlines!(ax, [radius]; color = :black, linestyle = :dash, label = "Radius √10 from complex poles")
axislegend(ax; position = :lt, labelsize = 11)
convergence_fig
