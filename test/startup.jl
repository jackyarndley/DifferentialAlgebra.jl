using Test, DifferentialAlgebra, LinearAlgebra
const DA = DifferentialAlgebra

# This is first in the fresh Pkg.test process, before any extension imports.
@testset "Fresh-process loading and optional diagonal fast path" begin
    @test !DA.isinitialized()
    @test Base.get_extension(DA, :DifferentialAlgebraIntervalArithmeticExt) === nothing
    @test Base.get_extension(DA, :DifferentialAlgebraOrdinaryDiffEqExt) === nothing
    @test_throws ArgumentError TaylorMethod(12)
    @test_throws MethodError TaylorMethod(12.0)
    @test !any(nameof(m) == :GenericLinearAlgebra for m in values(Base.loaded_modules))
    for T in (BigFloat, Rational{BigInt})
        x, = variables(T, 1; order = 2)
        values, vectors = DA.eigenpairs(Matrix(Diagonal([4 + x, 2 + 2x, 2 + 3x])))
        @test constant_term.(values) == [2, 2, 4]
        @test coefficient.(values, Ref([1])) == [2, 3, 1]
        @test constant_term.(vectors) == [0 0 1; 1 0 0; 0 1 0]
    end
end
