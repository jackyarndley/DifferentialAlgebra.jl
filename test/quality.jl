using Test, DifferentialAlgebra, Aqua

Aqua.test_all(DifferentialAlgebra; unbound_args = false)

@testset "Unbound arguments" begin
    # Julia 1.10 flags these generated outer constructors: an empty Vararg tuple
    # cannot infer its element type. Supported snapshots require nonempty domains
    # and patches. Exempt exactly these five methods, not their types or any newly
    # added overloads; newer Julia detectors already omit them.
    constructors = (
        DifferentialAlgebra.ModelCoordinates, PiecewisePolygonMap,
        ContinuousTaylorMap, TaylorModelPatch, PiecewiseTaylorModel,
    )
    exceptions = Set(only(methods(T)) for T in constructors)
    unbound = Test.detect_unbound_args(DifferentialAlgebra; recursive = true)
    @test isempty(setdiff(unbound, exceptions))
end
