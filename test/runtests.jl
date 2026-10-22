using SafeTestsets: @safetestset

include("support/fixtures.jl")
# Modules isolate names; all context-mutating cases remain sequential.
if ARGS == ["quality"]
    include("quality.jl")
else
    isempty(ARGS) || error("Expected no arguments or 'quality'")
    @safetestset "Startup" begin
        include("startup.jl")
    end
    @safetestset "Polynomial algebra" begin
        include("core.jl")
    end
    @safetestset "Analytic functions" begin
        include("functions.jl")
    end
    @safetestset "Maps and calculus" begin
        include("maps.jl")
    end
    @safetestset "Special functions" begin
        include("special.jl")
    end
    @safetestset "Contexts" begin
        include("context.jl")
    end
    @safetestset "Interval models" begin
        include("models.jl")
    end
    @safetestset "Adaptive maps" begin
        include("adaptive.jl")
    end
    @safetestset "Time evolution" begin
        include("evolution.jl")
    end
    @safetestset "Optional integrations" begin
        include("interop.jl")
    end

end
