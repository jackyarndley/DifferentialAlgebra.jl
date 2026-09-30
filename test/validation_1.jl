@testset verbose = true "TEST 1: DifferentialAlgebra initialization" begin
    @testset "1.1 Maximum order, number of variables, and number of monomials" begin
        for order = 1:10
            for nvar = 1:6
                @testset "order=$(order), nvar=$(nvar)" begin
                    DifferentialAlgebra.init(order, nvar)
                    @test order == DifferentialAlgebra.getMaxOrder()
                    @test nvar == DifferentialAlgebra.getMaxVariables()
                    @test factorial(order + nvar) / (factorial(order) * factorial(nvar)) == DifferentialAlgebra.getMaxMonomials()
                end
            end
        end
    end

    @testset "1.2 Truncation order" begin
        for k_user = 1:10
            @testset "1.2 Trunction order (k_user=$k_user)" begin
                DifferentialAlgebra.init(10, 6)

                DifferentialAlgebra.setTO(k_user)
                k_da = DifferentialAlgebra.getTO()
                @test k_da == k_user
            end
        end
    end

    @testset "1.3 Cutoff for the coefficients" begin
        for cutoff_user in (1e-14, 1e-15, 1e-16)
            @testset "1.3 Cutoff for the coefficients (cutoff_user=$cutoff_user)" begin
                DifferentialAlgebra.init(10, 6)
                DifferentialAlgebra.setEps(cutoff_user)
                cutoff_da = DifferentialAlgebra.getEps()
                @test cutoff_da == cutoff_user
            end
        end
    end
end
