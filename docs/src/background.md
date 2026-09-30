# Mathematical background

Differential algebra evaluates calculations in a truncated polynomial algebra
[Berz1999](@cite). For `n` independent perturbations and maximum total degree
`m`, a polynomial has the form

```math
p(\delta) = \sum_{|\alpha|\le m} c_\alpha\,\delta^\alpha,
\qquad |\alpha| = \alpha_1+\cdots+\alpha_n.
```

For a Taylor expansion of a smooth function about `a`,

```math
c_\alpha = \frac{1}{\alpha!}\,\partial^\alpha f(a).
```

Addition and multiplication propagate these coefficients, discarding terms above
the maximum degree. Analytic functions propagate their local series through
the same algebra. Differentiation and integration act on monomial exponents.

Replacing numeric initial conditions with polynomials lets a numerical
integrator propagate a local map of the solution. This is useful in nonlinear
uncertainty propagation [Valli2013](@cite); see the
[orbit integration example](generated/ode_integration.md). The integration error
and the error from truncating the polynomial are separate quantities.
