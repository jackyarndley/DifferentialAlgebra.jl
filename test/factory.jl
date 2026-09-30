@testset "Independent variables" begin
    x = variables(4; order = 4)
    @test x isa Vector{DA{Float64}}
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
