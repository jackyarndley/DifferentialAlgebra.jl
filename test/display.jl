@testset "Polynomial display" begin
    x, y = variables(2; order = 5)
    @test sprint(show, zero(x)) == "0.0"
    @test sprint(show, TaylorPolynomial(-2)) == "-2.0"
    @test sprint(show, x) == "1.0x₁"
    @test sprint(show, -y) == "-1.0x₂"
    @test sprint(show, x * y) == "1.0x₁x₂"
    @test sprint(show, 2 + 3x - y^2) == "2.0 + 3.0x₁ - 1.0x₂²"
    @test sprint(show, -2x + x * y - y^2) == "-2.0x₁ + 1.0x₁x₂ - 1.0x₂²"
    @test sprint(show, TaylorPolynomial(Inf)) == "Inf"
    @test sprint(show, TaylorPolynomial(-Inf)) == "-Inf"
    @test sprint(show, TaylorPolynomial(NaN)) == "NaN"
    @test sprint(show, TaylorPolynomial{Int}(typemin(Int))) == string(typemin(Int))
    @test sprint(show, TaylorPolynomial{Float32}(2)) == "2.0f0"
    @test sprint(show, Monomial(3.0, UInt32[1, 2])) == "Monomial(3.0, [1, 2])"
    @test sprint(show, MIME"text/plain"(), x + y) ==
        "TaylorPolynomial{Float64} polynomial in 2 variables (order ≤ 5):\n  1.0x₁ + 1.0x₂"
    @test occursin("x₁", sprint(show, MIME"text/plain"(), [x y; y x]))
    @test sprint(show, x + y; context = :compact => true) == "1.0x₁ + 1.0x₂"
    map = CompiledMap([x + y^2, y])
    @test sprint(show, map) == "CompiledMap{Float64}(2 outputs, 2 variables, degree 2, 4 nodes)"

    p = (1 + x + y)^5
    full = sprint(show, p)
    limited = sprint(show, p; context = :limit => true)
    compact = sprint(show, p; context = (:limit => true, :compact => true))
    @test !occursin('…', full) && occursin("x₂⁵", full)
    @test endswith(limited, " + …")
    @test endswith(compact, " + …")
    @test length(compact) < length(limited) < length(full)
    @test sprint(show, p) == full # Display never changes coefficient storage.

    variables(Rational{BigInt}, 1; order = 2)
    @test sprint(show, p) == "TaylorPolynomial(inactive)"
    @test sprint(show, MIME"text/plain"(), p) == "TaylorPolynomial(inactive)"
    @test sprint(show, variable(1, Rational{BigInt}) / 3) == "1//3x₁"
    setprecision(256) do
        c = BigFloat(1) / 3
        @test sprint(show, TaylorPolynomial{BigFloat}(c)) == sprint(show, c)
    end
end

@testset "Named variables and Unicode powers" begin
    x, y = variables((:x, :y); order = 12)
    p = 2 + 3x - x * y + y^12
    @test string(p) == "2.0 + 3.0x - 1.0xy + 1.0y¹²"
    @test coefficient(p, [0, 12]) == 1
    @test p([2, 1]) == 7
    @test string(differentiate(p, 2)) == "-1.0x + 12.0y¹¹"
    @test string(variable(2)) == "1.0y"
    @test occursin("xy", sprint(show, MIME"text/plain"(), [p]))

    q, r = variables(2; order = 12, names = ["q₁", "δr"])
    @test string(q^10 * r^2) == "1.0q₁¹⁰δr²"
    @test string(p) == "TaylorPolynomial(inactive)"
    @test string(copy(q + r)) == "1.0q₁ + 1.0δr"
    @test string(evaluate(CompiledMap([q + r^2]), [q, r])[1]) == "1.0q₁ + 1.0δr²"

    x = variables(12; order = 2)
    @test string(x[10] * x[12]) == "1.0x₁₀x₁₂"
    @test string(x[12]^2) == "1.0x₁₂²"
end
