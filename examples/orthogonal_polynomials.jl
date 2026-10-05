# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # An orthonormal Legendre basis
#
# Generate Pₙ by (n+1)Pₙ₊₁=(2n+1)xPₙ-nPₙ₋₁. Products in x and y
# form a total-degree basis on [-1,1]². Integrate each squared basis
# polynomial exactly in the truncated algebra, then normalize.
using DifferentialAlgebra
using LinearAlgebra
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

function legendre(x, degree)
    values = [one(x), x]
    for n in 1:(degree - 1)
        push!(values, ((2n + 1) * x * values[end] - n * values[end - 1]) / (n + 1))
    end
    return values[1:(degree + 1)]
end

function box_integral(p)
    for variable in 1:2
        primitive = integrate(p, variable)
        p = DifferentialAlgebra.substitute(primitive, variable, 1) -
            DifferentialAlgebra.substitute(primitive, variable, -1)
    end
    return constant_term(p)
end

# Products have degree six; leave one further order for antiderivatives.
x, y = variables((:x, :y); order = 7)
px, py = legendre(x, 3), legendre(y, 3)
indices = [(i, j) for i in 0:3 for j in 0:(3 - i)]
basis = [px[i + 1] * py[j + 1] for (i, j) in indices]
basis = [p / sqrt(box_integral(p^2)) for p in basis]
gram = [box_integral(p * q) for p in basis, q in basis]
error = maximum(abs, gram - I)
@assert error < 1.0e-13
println("Normalized total-degree basis: ", basis)
println("Basis Jacobian: ", jacobian(basis))
println("Maximum Gram matrix error: ", error)

#-
fig = Figure(size = (580, 470), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "Basis index", ylabel = "Basis index", title = "Integrated Gram matrix", aspect = DataAspect())
heat = heatmap!(ax, gram; colormap = :viridis, colorrange = (0, 1))
Colorbar(fig[1, 2], heat; label = "Inner product")
fig

# ## What the basis functions look like
# The signed surfaces show how mixed degrees resolve different shapes, while
# the Gram matrix above checks their orthonormality under the box integral.
grid = range(-1, 1; length = 65)
selected = (2, 3, 5, 6, 7, 9)
surfaces = [[basis[k]([a, b]) for a in grid, b in grid] for k in selected]
limit = maximum(v -> maximum(abs, v), surfaces)
basis_fig = Figure(size = (1110, 680), fontsize = 14)
for (panel, (k, values)) in enumerate(zip(selected, surfaces))
    i, j = indices[k]
    row, column = divrem(panel - 1, 3) .+ 1
    ax = Axis(basis_fig[row, column]; xlabel = "x", ylabel = "y", title = "Normalized L$i(x) L$j(y)", aspect = DataAspect())
    heatmap!(ax, grid, grid, values; colormap = :balance, colorrange = (-limit, limit))
end
Colorbar(basis_fig[1:2, 4]; colormap = :balance, limits = (-limit, limit), label = "Basis value (shared scale)")
basis_fig
