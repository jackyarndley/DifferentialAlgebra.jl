using Test, DifferentialAlgebra
using LinearAlgebra, GenericLinearAlgebra, ForwardDiff

import ForwardDiff: Dual

@testset verbose = true "Eigendecomposition" begin

    DifferentialAlgebra.initialize!(1, 6)

    A = [4.0 2.0 1.0; 2.0 3.0 1.0; 1.0 1.0 2.0]
    Ada = Matrix{TaylorPolynomial}(undef, 3, 3)
    Afd = Matrix{Dual{Nothing, Float64, 6}}(undef, 3, 3)

    k = 1
    p = zeros(6)
    for i in 1:3
        for j in 1:i
            p[k] = 1.0
            Ada[i, j] = A[i, j] + TaylorPolynomial(k, 1.0)
            Afd[i, j] = Dual(A[i, j], p...)
            p[k] = 0.0
            k += 1
        end
    end

    for i in 1:3
        for j in 1:(i - 1)
            Ada[j, i] = Ada[i, j]
            Afd[j, i] = Afd[i, j]
        end
    end

    Fda = LinearAlgebra.eigen(Hermitian(Ada))
    Ffd = LinearAlgebra.eigen(Hermitian(Afd))

    λda, Vda = Fda.values, Fda.vectors          # GenericLinearAlgebra + TaylorPolynomial
    λcc, Vcc = DifferentialAlgebra.eigenpairs(Ada)                   # Taylor eigenpairs in Julia
    λfd, Vfd = Ffd.values, Ffd.vectors          # GenericLinearAlgebra + ForwardDiff

    cλda = DifferentialAlgebra.constant_term.(λda)
    cλcc = DifferentialAlgebra.constant_term.(λcc)
    cλfd = [λ.value for λ in λfd]

    @testset verbose = true "Eigenvalues constant part" begin
        @test all(isapprox.(cλda, cλcc, atol = 1.0e-14, rtol = 1.0e-14))
        @test all(isapprox.(cλda, cλfd, atol = 1.0e-14, rtol = 1.0e-14))
    end

    Jλda = vcat([DifferentialAlgebra.linear_part.(λda)...]'...)
    Jλcc = vcat([DifferentialAlgebra.linear_part.(λcc)...]'...)
    Jλfd = vcat([[λ.partials...] for λ in λfd]'...)

    @testset verbose = true "Eigenvalues partials" begin
        @test all(isapprox.(Jλda, Jλcc, atol = 1.0e-13, rtol = 1.0e-13))
        @test all(isapprox.(Jλda, Jλfd, atol = 1.0e-14, rtol = 1.0e-14))
    end

    cVda = DifferentialAlgebra.constant_term.(Vda)
    cVcc = DifferentialAlgebra.constant_term.(Vcc)
    cVfd = [V.value for V in Vfd]

    scc = [sign.(cVda[1, i]) == sign.(cVcc[1, i]) ? 1 : -1 for i in 1:3]'
    sfd = [sign.(cVda[1, i]) == sign.(cVfd[1, i]) ? 1 : -1 for i in 1:3]'

    @testset verbose = true "Eigenvectors constant part" begin
        @test all(isapprox.(cVda, cVcc .* scc, atol = 1.0e-14, rtol = 1.0e-14))
        @test all(isapprox.(cVda, cVfd .* sfd, atol = 1.0e-14, rtol = 1.0e-14))
    end

    JVda = [[DifferentialAlgebra.linear_part(V)...] for V in Vda]
    JVcc = [[DifferentialAlgebra.linear_part(V)...] for V in Vcc]
    JVfd = [[V.partials...] for V in Vfd]

    @testset verbose = true "Eigenvectors partials" begin
        for i in 1:3
            for j in 1:3
                @test all(isapprox.(JVda[i, j], JVcc[i, j] .* scc[j], atol = 1.0e-13, rtol = 1.0e-13))
                @test all(isapprox.(JVda[i, j], JVfd[i, j] .* sfd[j], atol = 1.0e-13, rtol = 1.0e-13))
            end
        end
    end
end
