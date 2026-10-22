using Test, DifferentialAlgebra, SpecialFunctions
include("support/special_references.jl")

@testset "Independent Bessel coefficients at problematic centers" begin
    for (family, center, n, degrees, strings) in SPECIAL_REFERENCES
        # Each scalar type gets a distinct numerical case, rather than the
        # Cartesian product of every type, center, precision, degree and table.
        T = family == :K ? (center == 5 ? BigFloat : Float32) : Float64
        f = family == :I ? besselix : family == :K ? besselkx : besselj
        selected = T === BigFloat ? (1:2) : eachindex(degrees)
        for budget in (family === :I ? (0, 32 * 1024^2) : (0,))
            setprecision(128) do
                x, = variables(
                    T, 1; order = maximum(degrees[selected]), table_bytes = budget
                )
                result = f(n, T(center) + x)
                expected = parse.(BigFloat, strings)
                @test coefficient_type(result) === T
                @test all(
                    isapprox(
                        coefficient(result, [k]), T(c);
                        rtol = 1000eps(T)
                    ) for (k, c) in zip(degrees[selected], expected[selected])
                )
                if family == :I
                    negative = f(n, -T(center) + x)
                    @test coefficient(negative, [1]) ≈ -T(expected[2]) rtol = 1000eps(T)
                end
            end
        end
    end
end

@testset "BigFloat special functions and real singularities" begin
    setprecision(256) do
        x, = variables(BigFloat, 1; order = 3)
        # Independent 100-digit mpmath values, at exactly representable 17/8.
        i = parse(
            BigFloat,
            "0.809196566007602118254159830538643220568413967382672328252350891607" *
                "223341316661635156423216783493928"
        )
        k = parse(
            BigFloat,
            "0.209609303999835197149715854215400089359609810494235110566159288687" *
                "6975266704431100028063620394724142"
        )
        @test constant_term(besseli(2, big(17) / 8 + x)) ≈ i rtol = 1000eps(BigFloat)
        @test constant_term(besselk(2, big(17) / 8 + x)) ≈ k rtol = 1000eps(BigFloat)
        @test constant_term(polygamma(3, 2 + zero(x))) ≈
            6(zeta(BigFloat(4)) - 1) rtol = 100eps(BigFloat)
    end
    x, = variables(1; order = 2)
    for f in (bessely, besselk, besselkx)
        @test_throws DifferentialAlgebra.TaylorError f(2, -1 + x)
    end
    @test_throws ArgumentError polygamma(-1, x)
end
