module DifferentialAlgebra

import Random
using Printf: @printf
using LinearAlgebra: LinearAlgebra, I, lu, SingularException
import SpecialFunctions
import QuadGK

export DA, AlgebraicVector, AlgebraicMatrix, compiledDA, Monomial
export variables, variable, coefficient, constant_term, differentiate, integrate
export gradient, jacobian, hessian, evaluate, evaluate!, compile, CompiledMap, invert

include("engine.jl")
include("arithmetic.jl")
include("special_scalars.jl")
include("functions.jl")
include("coefficients.jl")
include("substitution.jl")
include("io.jl")
include("arrays.jl")
include("evaluation.jl")
include("linear_algebra.jl")
include("statistics.jl")
include("api.jl")
include("docs.jl")

end
