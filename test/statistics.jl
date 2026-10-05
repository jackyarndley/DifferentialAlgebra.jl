using Test, LinearAlgebra, DifferentialAlgebra

# Independent Cartesian enumeration, rather than the package's recursive basis
# traversal. Keep these small cases bounded: this deliberately enumerates a cube.
function expected_multiindices(order, dimension)
    indices = [
        collect(alpha) for alpha in Iterators.product(ntuple(_ -> 0:order, dimension)...)
            if sum(alpha) <= order
    ]
    return sort!(indices; by = alpha -> (sum(alpha), Tuple(-a for a in alpha)))
end

# Gaussian integration by parts (equivalently, Wick/Isserlis contractions):
# E[X_i X^beta] = mu_i E[X^beta] + sum_j beta_j Sigma[i,j] E[X^(beta-e_j)].
# Evaluate with exact rational means/covariances. This oracle never constructs a
# Taylor polynomial or differentiates a moment-generating function.
function gaussian_moment(alpha, mean, covariance)
    i = findfirst(!iszero, alpha)
    isnothing(i) && return one(eltype(mean))
    beta = collect(alpha)
    beta[i] -= 1
    value = mean[i] * gaussian_moment(beta, mean, covariance)
    for j in eachindex(beta)
        count = beta[j]
        iszero(count) && continue
        beta[j] -= 1
        value += count * covariance[i, j] * gaussian_moment(beta, mean, covariance)
        beta[j] += 1
    end
    return value
end

@testset verbose = true "Moments" begin
    @testset "Multi-indices" begin
        for (order, dimension) in ((0, 1), (3, 2), (6, 4), (3, 6))
            indices = DifferentialAlgebra.multiindices(order, dimension)
            expected = expected_multiindices(order, dimension)
            @test length(indices) == length(expected) == binomial(order + dimension, dimension)
            for (actual, alpha) in zip(indices, expected)
                @test actual == alpha
            end
        end
    end

    # Small inputs are inline; all 210 moments of each four-variable case are
    # computed from their definition in Julia instead of loaded from fixtures.
    cases = (
        (
            "Integer parameters",
            [3, 2, -1, 2],
            [4 2 1 1; 2 3 1 2; 1 1 2 1; 1 2 1 3],
            (1.0e-15, 1.0e-14),
        ),
        (
            "Floating parameters",
            [3.861346813186461, -1.487123146365428, -2.364168486116943, 1.0610402543840649],
            [
                4.564531354848048e-2 -2.003168404365469e-4 1.568048483265416e-3 -1.008404804051421e-3
                -2.003168404365469e-4 5.204846143436811e-1 1.003843436114855e-5 -1.008403740486084e-3
                1.568048483265416e-3 1.003843436114855e-5 3.2570408608486097e-2 2.10863483486048e-6
                -1.008404804051421e-3 -1.008403740486084e-3 2.10863483486048e-6 2.903406840638489e-1
            ],
            (1.0e-15, 5.0e-12),
        ),
    )

    @testset "Gaussian oracle identities" begin
        # Anchor the oracle itself with elementary exact moments and correlated
        # Wick contractions, including repeated indices and signed means.
        mean = Rational{BigInt}.([3, -2])
        covariance = Rational{BigInt}.([4 1; 1 9])
        @test gaussian_moment([0, 0], mean, covariance) == 1
        @test gaussian_moment([1, 0], mean, covariance) == 3
        @test gaussian_moment([2, 0], mean, covariance) == 13
        @test gaussian_moment([1, 1], mean, covariance) == -5
        @test gaussian_moment([4, 0], zero(mean), covariance) == 3 * 4^2
        @test gaussian_moment([2, 2], zero(mean), covariance) == 4 * 9 + 2 * 1^2
        @test gaussian_moment([3, 1], zero(mean), covariance) == 3 * 4 * 1
        @test gaussian_moment([6, 0], zero(mean), covariance) == 15 * 4^3
        @test gaussian_moment([1, 2], zero(mean), covariance) == 0
    end

    for (name, input_mean, input_covariance, tolerances) in cases
        @testset "$(name)" begin
            order, dimension = 6, 4
            mean = Rational{BigInt}.(input_mean)
            covariance = Rational{BigInt}.(input_covariance)
            expected_indices = expected_multiindices(order, dimension)

            @testset "Raw moments, centered=$(centered)" for centered in (true, false)
                mu = centered ? zero(mean) : mean
                x = variables(dimension; order)
                mgf = exp(Float64.(mu)' * x + x' * Float64.(covariance) * x / 2)
                indices, moments = DifferentialAlgebra.raw_moments(mgf, order)
                @test indices == expected_indices
                for (alpha, actual) in zip(indices, moments)
                    expected = Float64(gaussian_moment(alpha, mu, covariance))
                    @test actual ≈ expected atol = tolerances[centered ? 1 : 2] rtol = 1.0e-15
                end
            end

            @testset "Central moments" begin
                x = variables(dimension; order)
                mgf = exp(Float64.(mean)' * x + x' * Float64.(covariance) * x / 2)
                indices, moments = DifferentialAlgebra.central_moments(mgf, order)
                @test indices == expected_indices
                for (alpha, actual) in zip(indices, moments)
                    expected = Float64(gaussian_moment(alpha, zero(mean), covariance))
                    @test actual ≈ expected atol = 1.0e-11 rtol = 1.0e-15
                end
            end
        end
    end
end
