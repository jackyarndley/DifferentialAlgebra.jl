@testset "Independent variables" begin
    x = variables(4; order = 4)
    @test x isa Vector{TaylorPolynomial{Float64}}
    @test constant_term(x) == zeros(4)
    @test constant_term(jacobian(x)) == Matrix{Float64}(I, 4, 4)

    # Retrieve selected variables without invalidating an existing calculation.
    p = x[1] * x[3]
    selected = variable.([3, 4, 1, 1, 3])
    @test all(iszero, selected - x[[3, 4, 1, 1, 3]])
    @test p([2, 3, 4, 5]) == 8
    @test_throws ArgumentError variable(0)
    @test_throws ArgumentError variable(5)
end

@testset "Variable names" begin
    names = ["position", "velocity"]
    x, v = variables(Float32, names; order = 3)
    @test coefficient_type(x) === Float32
    @test string(x + v) == "1.0f0 position + 1.0f0 velocity"
    names[1] = "changed"
    @test string(x) == "1.0f0 position"

    # Invalid names must leave the existing algebra usable.
    @test_throws DimensionMismatch variables(2; order = 3, names = (:x,))
    @test_throws ArgumentError variables((:x, "x"); order = 3)
    @test_throws ArgumentError variables(("", "y"); order = 3)
    @test_throws ArgumentError variables(("x+y", "y"); order = 3)
    @test_throws ArgumentError variables((:x, 2); order = 3)
    @test_throws ArgumentError variables(2; order = 3, names = "xy")
    @test (x + v)([2, 3]) == 5

    initialize!(3, 2; names = (:r, :v))
    @test string(variable(1) + variable(2)) == "1.0 r + 1.0 v"
    initialize!(3, 2)
    @test string(variable(1) + variable(2)) == "1.0 x₁ + 1.0 x₂"
end
