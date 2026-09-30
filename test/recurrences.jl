@testset "Coefficient recurrences and composition agree" begin
    # A zero lookup budget selects Horner composition through polynomial products,
    # providing an independent algorithm for dense, sparse and short inputs.
    for T in (Float32, Float64, BigFloat)
        results = []
        for budget in (0, 32 * 1024^2)
            x, y, z = variables(T, 3; order = 6, table_bytes = budget)
            inputs = (
                one(x), 1 + x / 10, 2 + y^2 / 20, 2 + z^5 / 30,
                2 + x / 10 - y * z / 20, (1 + x / 10 + y / 20 - z / 30)^6,
            )
            functions = (
                exp, log, sin, cos, sinh, cosh, sqrt, inv,
                p -> p^T(0.3), p -> (1 + x - z^2) / p,
            )
            values = []
            for order in (3, 6), p in inputs, f in functions
                result = with_order(() -> f(p), order)
                push!(values, [coefficient(result, powers) for powers in DifferentialAlgebra.multiindices(6, 3)])
            end
            push!(results, values)
        end
        tolerance = T === Float32 ? T(2.0e-5) : T === Float64 ? T(2.0e-13) : T(2)^(-200)
        for (a, b) in zip(results...)
            @test maximum(abs, a - b) <= tolerance * max(one(T), maximum(abs, b))
        end
    end
    x, y = variables(Rational{BigInt}, 2; order = 8)
    h = x / 3 - y^2 / 7
    geometric = sum((-h)^k for k in 0:8)
    @test coefficient_norm(inv(1 + h) - geometric) == 0
    @test coefficient_norm((2 + y) / (1 + h) - (2 + y) * geometric) == 0
end
