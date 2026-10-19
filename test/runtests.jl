using SafeTestsets: @safetestset

# Modules isolate names, not the global algebra. Keep these tests sequential.
@safetestset "loading" begin
    include("loading.jl")
end
@safetestset "quality" begin
    include("quality.jl")
end
@safetestset "configuration" begin
    include("configuration.jl")
end
@safetestset "constructors" begin
    include("constructors.jl")
end
@safetestset "public api" begin
    include("public_api.jl")
end
@safetestset "api validation" begin
    include("api_validation.jl")
end
@safetestset "variables" begin
    include("variables.jl")
end
@safetestset "polynomial" begin
    include("polynomial.jl")
end
@safetestset "comparison operators" begin
    include("comparison_operators.jl")
end
@safetestset "analytic series" begin
    include("analytic_series.jl")
end
@safetestset "special functions" begin
    include("special_functions.jl")
end
@safetestset "coefficient types" begin
    include("coefficient_types.jl")
end
@safetestset "recurrences" begin
    include("recurrences.jl")
end
@safetestset "regressions" begin
    include("regressions.jl")
end
@safetestset "numerical edges" begin
    include("numerical_edges.jl")
end
@safetestset "linear algebra" begin
    include("linear_algebra.jl")
end
@safetestset "hessians" begin
    include("hessians.jl")
end
@safetestset "statistics" begin
    include("statistics.jl")
end
@safetestset "arrays" begin
    include("arrays.jl")
end
@safetestset "display" begin
    include("display.jl")
end
@safetestset "time series" begin
    include("time_series.jl")
end
@safetestset "sciml" begin
    include("sciml.jl")
end
@safetestset "taylor solver" begin
    include("taylor_solver.jl")
end
@safetestset "reused calculations" begin
    include("reused_calculations.jl")
end
@safetestset "domain splitting" begin
    include("domain_splitting.jl")
end
@safetestset "polygon ads" begin
    include("polygon_ads.jl")
end
@safetestset "interval models" begin
    include("interval_models.jl")
end
@safetestset "validated ads" begin
    include("validated_ads.jl")
end
@safetestset "interval polygon ads" begin
    include("interval_polygon_ads.jl")
end
@safetestset "continuous ads" begin
    include("continuous_ads.jl")
end
@safetestset "ads contracts" begin
    include("ads_contracts.jl")
end
@safetestset "tutorial tests" begin
    include("tutorial_tests.jl")
end
