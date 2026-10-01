# Moving from DACEjl

[DACEjl](https://github.com/arma1978/DACEjl) and DifferentialAlgebra.jl both
represent multivariate Taylor expansions truncated by total degree. Their APIs
and storage models differ. The [tutorial collection](../examples.md) provides
executable adaptations of DACEjl's 23 numbered examples.

| Operation in DACEjl | DifferentialAlgebra.jl |
|:--|:--|
| `init(order, n)`, then `DA(i)` | `x = variables(n; order)`, then `x[i]` |
| Named independent variables | `variables((:x, :y); order = 5)` |
| `DA(1.0)` | `TaylorPolynomial(1.0)` |
| `constant(p)` | `constant_term(p)` |
| `getcoefficient(p, powers)` | `coefficient(p, powers)` |
| `setcoefficient!(p, powers, value)` | `set_coefficient!(p, powers, value)` |
| `getmonomials(p)` | `monomials(p)` |
| `differentiate(p, i)`, `integrate(p, i)` | Same operations; arrays preserve their shape |
| `evaluate(p, point)` | `p(point)` or `evaluate(p, point)` |
| `compose(map, substitutions)` | `evaluate(map, substitutions)` |
| `invert(map)` | `invert(map)`; a nonsingular linear part is required |
| `DAVector`, `DAMatrix` | Ordinary `Vector` and `Matrix` with polynomial entries |
| Repeated map evaluation | `compiled = compile(map)`, then `compiled(point)` or `evaluate!` |
| Coefficient norms | `coefficient_norm(p)` gives the maximum absolute coefficient |
| Variable substitution | `DifferentialAlgebra.substitute(p, i, value)` |
| Variable translation/scaling | `DifferentialAlgebra.translate_variable(p, i, scale, shift)` |
| Working order | `with_order(() -> calculation(), order)` |
| Domain splitting | `adaptive_map(f, lower, upper; order, atol)` |
| Embedded ODE integrators | OrdinaryDiffEq solvers operating on polynomial states |

## Algorithms and numerical checks

Arithmetic uses a graded monomial basis with bounded lookup-table storage.
Elementary functions use coefficient recurrences or univariate series
composition. Inversion lifts a nonsingular linear inverse degree by degree.
Compiled maps share intermediate monomials between outputs. See
[polynomial maps](maps.md) and the [mathematical background](../background.md).

ADS recomputes the callback on each child domain. It combines discarded guard
orders, tail extrapolation and point checks to choose patches. Its estimates
are heuristic; validate independent points and inspect `map.converged`.
The [domain-splitting manual](domain-splitting.md) describes the acceptance rule,
resource limits and coordinate conventions. An ODE flow can be the callback;
there is no separate integrator-specific manifold type.

`DifferentialAlgebra.bounds(p)` estimates the polynomial range on the normalized
box [-1, 1]ⁿ by bounding each monomial. For another box, first compose with
`center + radius .* x`. This does not bound the truncation error of the original
function and does not use outward-rounded interval arithmetic.

Coefficient storage is parameterized by a concrete subtype of `Real`.
The operations used by a calculation must exist for that coefficient type;
transcendental functions may require floating-point operations. See
[coefficient types](coefficient-types.md) for precision and promotion rules.
Initialization defines a shared algebra, so initialize before concurrent
calculations; reinitialization invalidates existing uncompiled polynomials.

The examples check identities, inverse compositions, implicit-equation
residuals, orthogonality, independent ODE solutions and sampled ADS errors.
An elementary-function check at only the expansion center is insufficient:
it verifies the constant coefficient but says nothing about higher derivatives.
Likewise, agreement between ODE solvers should include polynomial coefficients
as well as nominal trajectories.
