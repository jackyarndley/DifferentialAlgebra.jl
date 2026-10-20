using Test, DifferentialAlgebra

# Deliberately unsupported axes, independent of OffsetArrays.
struct ShiftedVector{T} <: AbstractVector{T}
    data::Vector{T}
end

Base.size(v::ShiftedVector) = size(v.data)
Base.axes(v::ShiftedVector) = (0:(length(v) - 1),)
Base.getindex(v::ShiftedVector, i::Int) = v.data[i + 1]

@testset "Checked dimensions and coefficient storage" begin
    x, y = variables(2; order = 4)
    for counts in ([1], [1, 0, 0])
        for f in (coefficient, differentiate, integrate)
            @test_throws DimensionMismatch f(x, counts)
        end
        @test_throws DimensionMismatch set_coefficient!(copy(x), counts, 2)
        @test_throws DimensionMismatch monomial(counts)
    end
    for counts in ([-1, 0], [big(typemax(Int)) + 1, 0])
        for f in (coefficient, differentiate, integrate)
            @test_throws ArgumentError f(x, counts)
        end
        @test_throws ArgumentError set_coefficient!(copy(x), counts, 2)
    end
    counts = ShiftedVector([1, 0])
    for f in (coefficient, differentiate, integrate)
        @test_throws ArgumentError f(x, counts)
    end
    @test_throws ArgumentError set_coefficient!(copy(x), counts, 2)
    @test_throws ArgumentError Monomial(1, counts)
    @test_throws ArgumentError evaluate(x, ShiftedVector([0.1, 0.2]))
    for T in (Real, AbstractFloat)
        @test_throws ArgumentError TaylorPolynomial{T}()
        @test_throws ArgumentError variable(1, T)
        @test_throws ArgumentError TaylorPolynomial{T}(1)
        @test_throws ArgumentError TaylorPolynomial{T}(1, 1)
        @test_throws ArgumentError TaylorPolynomial{T}(T[1], 1, x.algebra)
        @test_throws ArgumentError TimeSeries(T[1, 2])
    end
    for T in (Int, Rational{BigInt}, Float32, Float64, BigFloat)
        p = TaylorPolynomial{T}(2) + TaylorPolynomial{T}(1, 3)
        q = float(p)
        @test typeof(q) === float(typeof(p))
        @test coefficient(q, [1, 0]) == 3
    end
    # Evaluation keeps its documented coordinate padding and ignored extras.
    p = 2 + x + 3y
    @test p([4.0]) == 6
    @test p([4.0, 5.0, 99.0]) == 21
    data = [1.0, 2.0, 3.0]
    s = TimeSeries(view(data, 1:2))
    data[1] = 9
    @test s.coefficients isa Vector{Float64}
    @test coefficient(s, 0) == 1
    @test TimeSeries(1:3).coefficients == [1, 2, 3]
    for f in (round, trunc), p in (x + 2.3, TaylorPolynomial(Inf))
        @test_throws MethodError f(p)
        @test_throws MethodError f(Int, p)
    end
end

@testset "Contextual calculus and map errors" begin
    x, y = variables(2; order = 3)
    @test_throws r"Variable index 3 is outside 1:2" differentiate(x, 3)
    @test_throws r"Variable index 0 is outside 1:2" integrate(y, 0)
    @test_throws r"Variable power -1 is outside 0:3" (
        DifferentialAlgebra.divide_variable(x, 1, -1)
    )
    map = CompiledMap([x + y^2, y])
    @test_throws r"Output has 1 entries; expected 2" (
        evaluate!(zeros(1), map, [0.1, 0.2], zeros(3))
    )
    @test_throws r"Workspace needs at least 3 entries" (
        evaluate!(zeros(2), map, [0.1, 0.2], zeros(2))
    )
    @test_throws r"singular linear part" invert([x, 2x])
end
