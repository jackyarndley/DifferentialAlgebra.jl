using Test, DifferentialAlgebra

import DiffEqBase, OrdinaryDiffEqCore

@testset "Tangent series and recorded integrals" begin
    s = TimeSeries([0.0, 1.0, zeros(6)...])
    for (f, signs) in ((tan, ones(8)), (tanh, [(-1)^((k - 1) ÷ 2) for k in 1:8]))
        expected = [0, 1, 0, 1 / 3, 0, 2 / 15, 0, 17 / 315] .* signs
        @test f(s).coefficients ≈ expected
        expansion = only(taylor_expand((u, p, t) -> [f(t)], [0.0], 0.0; order = 8))
        @test expansion.coefficients ≈ [0; expected ./ (1:8)]
    end
end

@testset "Direct dense Horner interpolation" begin
    ext = Base.get_extension(DifferentialAlgebra, :DifferentialAlgebraOrdinaryDiffEqExt)
    cache = ext.TaylorCache(3, nothing, nothing, zeros(2))
    # Independent polynomials: 1+2t+3t²+4t³ and 5-2t².
    k = [[1.0, 5.0], [2.0, 0.0], [3.0, -2.0], [4.0, 0.0]]
    initial = first(k)
    expected = ([1.978, 4.82], [4.88, -1.2], [13.2, -4], [24.0, 0.0], [0.0, 0.0])
    for D in 0:4, idxs in (nothing, [2, 1], 2, view([1, 2], 1:1))
        indices = idxs === nothing ? [1, 2] : idxs isa Integer ? [idxs] : idxs
        out = fill(NaN, length(indices))
        wanted = expected[D + 1][indices]
        OrdinaryDiffEqCore._ode_interpolant!(out, 0.5, 0.6, initial, initial, k, cache, idxs, Val{D}, nothing)
        @test out ≈ wanted
        actual = OrdinaryDiffEqCore._ode_interpolant(0.5, 0.6, initial, initial, k, cache, idxs, Val{D}, nothing)
        @test (idxs isa Integer ? [actual] : actual) ≈ wanted
    end
    @test_throws DimensionMismatch OrdinaryDiffEqCore._ode_interpolant!(zeros(1), 0.5, 0.6, initial, initial, k, cache, nothing, Val{0}, nothing)
    @test_throws ArgumentError OrdinaryDiffEqCore._ode_interpolant!(zeros(2), 0.5, 0.6, initial, initial, k, cache, nothing, Val{-1}, nothing)
end
