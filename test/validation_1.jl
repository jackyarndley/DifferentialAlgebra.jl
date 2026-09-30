@testset verbose = true "TEST 1: DifferentialAlgebra initialization" begin
    @testset "1.1 Maximum order, number of variables, and number of monomials" begin
        for order in 1:10
            for nvar in 1:6
                @testset "order=$(order), nvar=$(nvar)" begin
                    DifferentialAlgebra.initialize!(order, nvar)
                    @test order == DifferentialAlgebra.max_order()
                    @test nvar == DifferentialAlgebra.nvariables()
                    @test factorial(order + nvar) / (factorial(order) * factorial(nvar)) == DifferentialAlgebra.nmonomials()
                end
            end
        end
    end

    @testset "1.2 Truncation order" begin
        for k_user in 1:10
            @testset "1.2 Trunction order (k_user=$k_user)" begin
                DifferentialAlgebra.initialize!(10, 6)

                DifferentialAlgebra.set_truncation_order!(k_user)
                k_da = DifferentialAlgebra.truncation_order()
                @test k_da == k_user
            end
        end
    end

    @testset "1.3 Cutoff for the coefficients" begin
        for cutoff_user in (1.0e-14, 1.0e-15, 1.0e-16)
            @testset "1.3 Cutoff for the coefficients (cutoff_user=$cutoff_user)" begin
                DifferentialAlgebra.initialize!(10, 6)
                DifferentialAlgebra.set_coefficient_tolerance!(cutoff_user)
                cutoff_da = DifferentialAlgebra.coefficient_tolerance()
                @test cutoff_da == cutoff_user
            end
        end
    end
end
