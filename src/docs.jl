@doc """
    DA{T}(c = zero(T))
    DA(c::Real)
    DA(i::Integer, c::Real)

A multivariate Taylor polynomial with coefficients of concrete real type `T`.
`DA(c)` creates a constant with coefficient type `typeof(float(c))`.
`DA{T}(c)` explicitly chooses storage type `T`. The two-argument constructor
creates `c` times variable `i`; index zero creates a constant.

Use [`variables`](@ref) to initialize the algebra. A polynomial is callable:
`p(point)` evaluates it numerically or composes it with polynomial coordinates.
""" DA

@doc """
    variable(i, T = Float64)

Return independent variable `i` in the current algebra with coefficient type `T`.
Indices start at one. This does not reinitialize or invalidate other polynomials.
""" variable

@doc """
    init(order, n; table_bytes = 32 * 1024^2)

Initialize the global algebra with maximum total degree `order` and `n` variables.
Existing polynomials become invalid. Both dimensions must be positive.
Use [`variables`](@ref) to initialize and construct the variables in one call.
`table_bytes` bounds lookup tables, excluding polynomial storage and basis metadata.
Configure the algebra before launching concurrent calculations.
""" init

@doc """
    coefficient(p, exponents)

Return the coefficient of a monomial, without factorial scaling.
Supply exactly one nonnegative exponent per independent variable, with total
degree at most the initialized order. Throws on invalid dimensions or exponents.
""" coefficient

@doc """
    integrate(p, i)
    integrate(p, counts)

Integrate a polynomial with respect to variable `i`, choosing zero integration constant.
A vector `counts` requests repeated integrals. Terms above the working order are
discarded. The coefficient type follows the scalar division operation.
An array and an integer index integrate each entry, preserving the array's shape.
""" integrate

@doc """
    gradient(p)

Return the vector of polynomial partial derivatives of `p`.
Use `constant_term(gradient(p))` for the gradient at the expansion point.
""" gradient

@doc """
    jacobian(p)
    jacobian(map)

Return a matrix of polynomial partial derivatives.
Rows correspond to map components and columns to independent variables.
A scalar polynomial produces a one-row matrix.
""" jacobian

@doc """
    hessian(p)
    hessian(map)

Return the matrix of second polynomial partial derivatives of `p`.
For a vector of polynomials, return a vector of Hessian matrices.
""" hessian

@doc """
    CompiledMap(p)
    CompiledMap(polynomials)

Construct a shared evaluation tree for repeated evaluation or composition.
Call `map(point)` or [`evaluate`](@ref) to evaluate it. Numeric evaluation remains
valid after algebra reinitialization because the map owns its coefficients.
Polynomial composition requires its original algebra.

`getDim`, `getOrd`, `getVars` and `getTerms` inspect its output dimension,
maximum degree, highest variable index and number of evaluation-tree nodes.
""" CompiledMap

@doc """
    compile(p)
    compile(polynomials)

Build a [`CompiledMap`](@ref) from a polynomial or vector of polynomials.
""" compile

@doc """
    evaluate(p, coordinates)
    evaluate(map, coordinates)
    evaluate(map, coordinates, result)

Evaluate a polynomial, vector of polynomials or compiled map.
Numeric coordinates return values; polynomial coordinates compose polynomials.
Missing coordinates are zero and extra coordinates are ignored. A scalar
polynomial returns a scalar; a map returns a vector. The three-argument form
writes into `result`. `p(coordinates)` and `map(coordinates)` provide call syntax.
""" evaluate

@doc """
    evaluate!(result, map::CompiledMap, coordinates, work)

Evaluate a compiled map numerically into reusable buffers.
`result` must have `DifferentialAlgebra.getDim(map)` entries and `work` at least
`DifferentialAlgebra.getOrd(map)+1` entries. Output, workspace and inputs must
not overlap. The output and workspace have the same scalar element type.
Scalar arithmetic may allocate even when these buffers are reused.
Use [`evaluate`](@ref) for polynomial composition.
""" evaluate!

@doc """
    invert(map)

Compute the local inverse of a polynomial map by degree-by-degree lifting.
The linear part must be nonsingular. The returned polynomial includes the
shift by the map's constant value and accepts output coordinates.
A partial map is completed with identity coordinates before inversion.
""" invert

@doc """
    Monomial(coefficient, exponents)

A coefficient and its vector of nonnegative integer exponents.
The exponents are stored as `UInt32` values. Retrieve monomials from a polynomial
with `DifferentialAlgebra.getMonomials(p)`.
""" Monomial

@doc """
    getCoefficient(p, exponents)
    getCoefficient(m::Monomial)

Read a coefficient using permissive exponent dimensions.
Missing exponents are zero, extra exponents are ignored, and total degrees above
the initialized order return zero. Prefer [`coefficient`](@ref) for explicit
dimension validation. `setCoefficient!` changes a coefficient in place.
""" getCoefficient

@doc """
    norm(p::DA, type = 0)

Compute a norm over all polynomial coefficients.
Type zero returns the maximum absolute coefficient; type one returns their
absolute sum; higher integer types use the corresponding power norm.
Unlike scalar polynomial comparisons, this includes nonconstant coefficients.
""" norm

@doc """
    getRawMoments(mgf, order)

Return multi-indices and raw moments through total degree `order`.
`mgf` must be a moment-generating function about zero. Coefficients are multiplied
by the factorials of their exponents to recover the moments.
""" getRawMoments

@doc """
    getCentralMoments(mgf, order)

Return multi-indices and moments after centering the moment-generating function.
The input must have constant coefficient one.
""" getCentralMoments
