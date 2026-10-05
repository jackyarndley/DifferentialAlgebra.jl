# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Linear algebra with polynomial entries
#
# Factorization, solving, determinants and traces use Julia's LinearAlgebra.
# The constant matrix must be nonsingular for a local polynomial inverse.
using DifferentialAlgebra
using LinearAlgebra
using CairoMakie

x, y = variables((:x, :y); order = 3)
A = [1 + 0.2x x - y; 0.5y 2 - 0.1x + y]
B = [x one(x); y x + y]
b = [1 + x, -0.5 + y]
F = lu(A)
z = F \ b
Ainverse = F \ Matrix{typeof(x)}(I, 2, 2)
solve_error = maximum(coefficient_norm, A * z - b)
inverse_error = maximum(coefficient_norm, A * Ainverse - Matrix{typeof(x)}(I, 2, 2))
@assert max(solve_error, inverse_error) < 1.0e-13
println("Solution of A z = b: ", z)
println("trace(A B) = ", tr(A * B))
println("det(A) = ", det(A))
println("Frobenius norm of A B: ", sqrt(sum(abs2, A * B)))
println((solve_error = solve_error, inverse_error = inverse_error))

# Julia's ordinary arrays preserve matrix-product associativity in the
# truncated algebra. Numeric evaluation need not commute with multiplication:
# powers above the working order have already been discarded.
associativity_error = maximum(coefficient_norm, (A * B) * b - A * (B * b))
@assert associativity_error < 1.0e-13
println("Matrix-product associativity residual: ", associativity_error)

# Compose each matrix entry with a shifted coordinate system.
shifted = [evaluate(p, [x + 0.1, y - 0.2]) for p in A]
println("A(x+0.1, y-0.2) = ", shifted)

# ## Solving the stored matrix at physical points
# Polynomial LU solves coefficient equations through order three. Numeric
# evaluation of that truncated solution has a residual away from the center.
# Compare it with a fresh scalar LU solve; conditioning is a separate issue.
matrix(v) = [1 + 0.2v[1] v[1] - v[2]; 0.5v[2] 2 - 0.1v[1] + v[2]]
rhs(v) = [1 + v[1], -0.5 + v[2]]
compiled = compile(z)
xs, ys = range(-0.8, 0.8; length = 65), range(-0.6, 0.6; length = 65)
residuals = [norm(matrix([a, b]) * compiled([a, b]) - rhs([a, b]), Inf) for a in xs, b in ys]
conditioning = [cond(matrix([a, b])) for a in xs, b in ys]
fig = Figure(size = (1230, 420), fontsize = 14)
ax = Axis(fig[1, 1]; xlabel = "physical x (y=0)", ylabel = "solution component", title = "Scalar solve versus order-three map")
for j in 1:2
    lines!(ax, xs, [(matrix([a, 0]) \ rhs([a, 0]))[j] for a in xs]; label = "Scalar z$j", linewidth = 2)
    lines!(ax, xs, [compiled([a, 0])[j] for a in xs]; linestyle = :dash, label = "Taylor z$j")
end
axislegend(ax; position = :rt, labelsize = 11)
ax = Axis(fig[1, 2]; xlabel = "physical x", ylabel = "physical y", title = "log₁₀ numeric solve residual", aspect = DataAspect())
heat = heatmap!(ax, xs, ys, log10.(max.(residuals, eps(Float64))); colormap = :magma)
Colorbar(fig[1, 3], heat)
ax = Axis(fig[1, 4]; xlabel = "physical x", ylabel = "physical y", title = "Matrix condition number", aspect = DataAspect())
heat = heatmap!(ax, xs, ys, conditioning; colormap = :viridis)
Colorbar(fig[1, 5], heat)
fig

# These are sampled residuals/condition numbers. A nonsingular constant matrix
# establishes a local formal inverse, not a certified inverse throughout a box.
