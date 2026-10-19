using Test, DifferentialAlgebra, LinearAlgebra

@test !any(m -> nameof(m) == :GenericLinearAlgebra, values(Base.loaded_modules))
for T in (BigFloat, Rational{BigInt})
    x, = variables(T, 1; order = 2)
    diagonal = [4 + x, 2 + 2x, 2 + 3x]
    A = Matrix(Diagonal(diagonal))
    values, vectors = DifferentialAlgebra.eigenpairs(A)
    @test coefficient_type(eltype(values)) === BigFloat
    @test coefficient_type(eltype(vectors)) === BigFloat
    @test constant_term.(values) == BigFloat[2, 2, 4]
    @test coefficient.(values, Ref([1])) == BigFloat[2, 3, 1]
    @test constant_term.(vectors) == BigFloat[0 0 1; 1 0 0; 0 1 0]
    @test all(degree(v) == 0 for v in vectors)
    values[1].coeffs[1] = 7
    @test constant_term(A[2, 2]) == 2
    vectors[1, 1].coeffs[1] = 7
    @test constant_term(vectors[1, 2]) == 0
end
@test !any(m -> nameof(m) == :GenericLinearAlgebra, values(Base.loaded_modules))
