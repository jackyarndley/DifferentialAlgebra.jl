using Test, DifferentialAlgebra

@test !DifferentialAlgebra.isinitialized()
@test_throws ArgumentError max_order()
time_expansion = only(taylor_expand((u, p, t) -> u, [1.0], 0.0; order = 12))
@test time_expansion(0.1) ≈ exp(0.1)
@test !DifferentialAlgebra.isinitialized()
initial_map = adaptive_map(x -> 1 + x[1], [-1.0], [1.0])
@test !DifferentialAlgebra.isinitialized()
@test initial_map([0.3]) ≈ 1.3
