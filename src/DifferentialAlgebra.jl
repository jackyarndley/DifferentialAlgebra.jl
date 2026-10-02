module DifferentialAlgebra

import Random
using LinearAlgebra: LinearAlgebra, I, lu
import SpecialFunctions
import QuadGK

export TaylorPolynomial, Monomial
export variables, variable, coefficient, constant_term, differentiate, integrate
export gradient, jacobian, hessian, evaluate, evaluate!, compile, CompiledMap, invert
export coefficient_type, set_coefficient!, monomials, monomial, exponents, degree
export coefficient_norm, linear_part, nvariables, max_order, noutputs, nnodes
export initialize!, truncation_order, set_truncation_order!, with_order
export coefficient_tolerance, set_coefficient_tolerance!
export adaptive_map, adaptive_flow, PiecewiseTaylorMap, TaylorPatch
export GuardedTail, ExtrapolatedTail, LastTerms
export TimeSeries, taylor_expand, TaylorMethod

include("basis.jl")
include("polynomial.jl")
include("arithmetic.jl")
include("special_scalars.jl")
include("functions.jl")
include("coefficients.jl")
include("substitution.jl")
include("display.jl")
include("evaluation.jl")
include("linear_algebra.jl")
include("statistics.jl")
include("ads_estimators.jl")
include("domain_splitting.jl")
include("adaptive_flow.jl")
include("time_series.jl")
include("time_recurrence.jl")
include("precompile.jl")

end
