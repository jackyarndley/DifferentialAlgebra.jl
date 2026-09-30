@testset "Polynomial display" begin
    x, y = variables(2; order = 5)
    @test sprint(show, zero(x)) == "0.0"
    @test sprint(show, DA(-2)) == "-2.0"
    @test sprint(show, x) == "x1"
    @test sprint(show, -y) == "-x2"
    @test sprint(show, x * y) == "x1*x2"
    @test sprint(show, 2 + 3x - y^2) == "2.0 + 3.0*x1 - x2^2"
    @test sprint(show, -2x + x * y - y^2) == "-2.0*x1 + x1*x2 - x2^2"
    @test sprint(show, DA(Inf)) == "Inf"
    @test sprint(show, DA(-Inf)) == "-Inf"
    @test sprint(show, DA(NaN)) == "NaN"
    @test sprint(show, DA{Int}(typemin(Int))) == string(typemin(Int))
    @test sprint(show, DA{Float32}(2)) == "2.0f0"
    @test sprint(show, Monomial(3.0, UInt32[1, 2])) == "Monomial(3.0, [1, 2])"
    @test sprint(show, MIME"text/plain"(), x + y) ==
        "DA{Float64} polynomial in 2 variables (order ≤ 5):\n  x1 + x2"
    @test occursin("x1", sprint(show, MIME"text/plain"(), [x y; y x]))
    @test sprint(show, x + y; context = :compact => true) == "x1 + x2"

    p = (1 + x + y)^5
    full = sprint(show, p)
    limited = sprint(show, p; context = :limit => true)
    compact = sprint(show, p; context = (:limit => true, :compact => true))
    @test !occursin('…', full) && occursin("x2^5", full)
    @test endswith(limited, " + …")
    @test endswith(compact, " + …")
    @test length(compact) < length(limited) < length(full)
    @test sprint(show, p) == full # Display never changes coefficient storage.

    variables(Rational{BigInt}, 1; order = 2)
    @test sprint(show, p) == "DA(inactive)"
    @test sprint(show, MIME"text/plain"(), p) == "DA(inactive)"
    @test sprint(show, variable(1, Rational{BigInt}) / 3) == "1//3*x1"
    setprecision(256) do
        c = BigFloat(1) / 3
        @test sprint(show, DA{BigFloat}(c)) == sprint(show, c)
    end
end
